#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ectrans4py:

A Python interface to spectral transforms from ecTrans, using cTypesForFortran for the
Fortran/Python binding.

Two interfaces are provided, sharing one geometry inquiry (``trans_inq4py``) and
one set of Legendre assets (``get_legendre_assets``):

* Serial (single task): the LAM/Gaussian wrappers ``sp2gp_*``/``gp2sp_*``,
  ``etrans_inq4py``, ``get_legendre_assets``, and ``trans_inq4py`` with
  ``KRESOL<=0`` (it self-initialises from the grid parameters).
* Distributed-memory (MPI): initialise with ``mpl_init4py`` (attaches to an
  existing communicator, e.g. mpi4py), or with ``mpl_init_comm4py`` to run ecTrans on
  the tasks of a communicator of one's own (e.g. split from ``MPI_COMM_WORLD``) only,
  set up the processor grid/resolution with
  ``setup_trans0_4py``/``setup_trans_4py`` (``LDSPLIT=.TRUE.``), query the local
  geometry with ``trans_inq4py`` passing the returned ``KRESOL``, transform local
  arrays with the ``*_dist4py`` routines, move data with ``dist_spec``/
  ``gath_spec``/``dist_grid``/``gath_grid``, and compute global norms with
  ``specnorm4py``/``gpnorm_trans4py``. With several V-sets (``KPRTRW < NPROC``), the
  spectral fields are distributed over the V-sets as well: use the ``*_vset4py``
  bindings, with a map from each field to its V-set (``trans_inq_vset4py`` returns
  the task's). ``inv_trans_all_dist4py``/``dir_trans_all_dist4py`` transform winds
  and scalar fields together, in a single call. The repacking of the spectral arrays
  between the layout of the callers and that of ecTrans is blocked and threaded with
  OpenMP when the environment variable ``ECTRANS4PY_OMP`` is set (to ``1``, say), with
  bit-identical results (``spec_layout_inq4py``).

Several resolutions can coexist (``KMAX_RESOL`` of ``setup_trans0_4py``): every
distributed transform, scatter/gather and norm takes a trailing ``KRESOL``, the tag
``setup_trans_4py`` returned (default 1, the first resolution set up), and runs on
that resolution. A ``KRESOL`` that is not a resolution set up by ``setup_trans_4py``
raises ``RuntimeError`` with error code 101, instead of ecTrans aborting the job. The
bindings take positional arguments only (the ctypesForFortran wrappers forward
``*args``), so pass ``KRESOL`` positionally.

A resolution can be set up for grid-point space only (``setup_trans_4py`` with
``LDGRIDONLY`` True, CPU only), as the IFS does for its radiation grid: its grid-point
partition is distributed (``dist_grid4py``/``gath_grid4py``) and inquired
(``trans_inq_grid4py``, with the region counts of ``setup_trans0_regions4py``), and
every binding that needs spectral or Fourier space raises ``RuntimeError`` with error
code 102 for it, instead of ecTrans crashing.
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import os
import resource
import numpy as np
import ctypesForFortran
from ctypesForFortran import addReturnCode, treatReturnCode, IN, OUT, array2string
import platform

# Shared objects library
########################
system = platform.system()
if system == "Linux":
    platform_ext = "so"
elif system == "Darwin":
    platform_ext = "dylib"
else:
    raise NotImplementedError("ectrans4py does not support Windows")

# Local names of dp and sp libraries in the directory
lib_basenames = [f"libectrans4py_{p}.{platform_ext}" for p in ["dp", "sp"]]
LD_LIBRARY_PATH = [p for p in os.environ.get('LD_LIBRARY_PATH', '').split(':') if p != '']
lpath = LD_LIBRARY_PATH + [
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lib'),
    os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lib64'),
]


def _find_library(basename):
    for d in lpath:
        candidate = os.path.join(d, basename)
        if os.path.exists(candidate):
            return candidate
    return None


for _p in ["dp", "sp"]:
    shared_objects_library = _find_library(f"libectrans4py_{_p}.{platform_ext}")
    if shared_objects_library is not None:
        _prec = _p
        break
else:
    msg = (
        f'libectrans4py_{{dp/sp}}.{platform_ext} was not found in any of '
        f'potential locations: {str(lpath)}. '
        f'You can specify a different location using env var LD_LIBRARY_PATH'
    )
    raise FileNotFoundError(msg)

# Floating-point type of the field data (spectral / grid-point) crossing the
# Fortran interface (the JPRB arrays): float32 for the single-precision build,
# float64 for the double. So are the Gaussian weights, Legendre polynomials and
# resolution deltas; only the sines of the Gaussian latitudes (PMU), which
# ecTrans returns in double precision, are float64 regardless.
_REAL = np.float32 if _prec == "sp" else np.float64

ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(shared_objects_library)

# Initialization
################

def init_env(omp_num_threads=None,
             no_mpi=True,
             unlimited_stack=True,
             ):
    """
    Set adequate environment for the inner libraries.

    :param int omp_num_threads: sets OMP_NUM_THREADS
    :param bool no_mpi: environment variable DR_HOOK_NOT_MPI set to 1
    :param unlimited_stack: equivalent to 'ulimit -s unlimited'
    """
    # because arpifs library is compiled with MPI & openMP
    if omp_num_threads is not None:
        os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)
    if no_mpi:
        os.environ['DR_HOOK_NOT_MPI'] = '1'
    if unlimited_stack:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))

# Transforms interfaces
#######################

@array2string(0)
@ctypesFF()
def ectrans_version():
    """
    Return the version string of ecTrans.

    Returns:\n
    1) CD_VERSION_STRING: version string of ecTrans (always 14 elements so must be trimmed)
    """
    return ([], [(str, (1, 14), OUT)], None)

@treatReturnCode
@ctypesFF()
@addReturnCode
def get_legendre_assets(KSIZEJ, KTRUNC, KSLOEN, KSPOLEGL, KLOEN, KNUMMAXRESOL):
    """
    Fetch arrays relevant for performing the Legendre transform.
    KNMENG and PGW are specified across the full globe, pole to pole. PRPNM is specified across the
    Northern hemisphere only.

    Args:\n
    1) KSIZEJ: number of latitudes in grid-point space
    2) KTRUNC: truncation
    3) KSLOEN: Size of KLOEN
    4) KSPOLEGL: the second dimension of the array storing all of the Legendre polynomials, equal to
       sum([truncation + 2 - im for im in range(truncation+1)])
    5) KLOEN: number of points on each latitude row
    6) KNUMMAXRESOL: maximum number of troncatures handled

    Returns:\n
    1) KNMENG: cut-off zonal wavenumber
    2) PGW: Gaussian weights
    3) PMU: sines of the Gaussian latitudes
    4) PRPNM: associated Legendre polynomials
    """
    return ([KSIZEJ, KTRUNC, KSLOEN, KSPOLEGL, KLOEN, KNUMMAXRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), OUT),
             (_REAL, (KSLOEN,), OUT),
             (np.float64, (KSLOEN,), OUT),
             (_REAL, (KSLOEN//2,KSPOLEGL), OUT)],
            None)

@treatReturnCode
@ctypesFF()
@addReturnCode
def etrans_inq4py(KSIZEI, KSIZEJ,
                 KPHYSICALSIZEI, KPHYSICALSIZEJ,
                 KTRUNCX, KTRUNCY,
                 KNUMMAXRESOL,
                 PDELATX, PDELATY):
    """
    Simplified wrapper to ETRANS_INQ.

    Args:\n
    1,2) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    3,4) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    5,6) KTRUNCX, KTRUNCY: troncatures
    7) KNUMMAXRESOL: maximum number of troncatures handled
    8,9) PDELTAX, PDELTAY: resolution along x,y axis

    Returns:\n
    1) KGPTOT: number of gridpoints
    2) KSPEC: number of spectral coefficients
    """
    return ([KSIZEI, KSIZEJ,
             KPHYSICALSIZEI, KPHYSICALSIZEJ,
             KTRUNCX, KTRUNCY,
             KNUMMAXRESOL,
             PDELATX, PDELATY],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, None, IN),
             (_REAL, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def trans_inq4py(KRESOL, KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL):
    """
    Wrapper to TRANS_INQ: extract the geometry of a resolution.

    Serves both the serial and the distributed setup through KRESOL:
    * KRESOL <= 0: self-initialise serially (single task) from the grid parameters
      (KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL); local sizes equal global.
    * KRESOL >  0: inquire the resolution already set up by setup_trans_4py (after the
      parallel setup_trans0_4py); the local sizes are this task's partition and the
      grid parameters are unused. A KRESOL that is not a resolution set up by
      setup_trans_4py raises RuntimeError with error code 101.

    Args:\n
    1) KRESOL: resolution tag from setup_trans_4py, or <=0 to self-initialise
    2) KSIZEJ: number of Gaussian latitudes (sizes KNMENG, PMU, PGW)
    3) KTRUNC: troncature                         (self-init only)
    4) KSLOEN: Size of KLOEN                       (self-init only)
    5) KLOEN: number of points on each latitude row (self-init only)
    6) KNUMMAXRESOL: maximum number of troncatures handled (self-init only)

    Returns:\n
    1) KGPTOT: local number of gridpoints
    2) KSPEC: local number of spectral coefficients
    3) KSPEC2: local number of (doubled) spectral coefficients
    4) KGPTOTG: global number of gridpoints
    5) KSPEC2G: global number of (doubled) spectral coefficients
    6) KSMAX: spectral truncation T
    7) KNMENG: cut-off zonal wavenumber (per latitude)
    8) PMU: sines of the Gaussian latitudes (global, length KSIZEJ), float64 whatever the
       precision of the build
    9) PGW: Gaussian weights (global, length KSIZEJ)
    """
    return ([KRESOL, KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, (KSIZEJ,), OUT),
             (np.float64, (KSIZEJ,), OUT),
             (_REAL, (KSIZEJ,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def trans_inq_grid4py(KRESOL, KDGL, KNRNS, KNREW):
    """
    Wrapper to TRANS_INQ for the grid-point partition tables of a distributed resolution,
    spectral or grid-point-only (setup_trans_4py LDGRIDONLY), as the IFS inquires them for
    its radiation grid.

    Args:\n
    1) KRESOL: resolution tag from setup_trans_4py
    2) KDGL: number of latitude rows of that resolution
    3,4) KNRNS, KNREW: numbers of N-S bands and of E-W regions of the grid-point partition,
       k_regions_ns and k_regions_ew of setup_trans0_4py (or setup_trans0_regions4py)

    Returns:\n
    1) NGPTOT: grid points on this task
    2) NGPTOTG: grid points on the globe
    3) NGPTOTMX: largest number of grid points on a task
    4) NFRSTLOFF: offset of the first latitude of this task's N-S band
    5) NPTRFLOFF: offset of the row of that latitude in NSTA/NONL
    6,7) MY_REGION_NS, MY_REGION_EW: this task's region (1-based)
    8) NFRSTLAT (KNRNS,): first latitude of each N-S band
    9) NLSTLAT (KNRNS,): last latitude of each N-S band
    10) NPTRFRSTLAT (KNRNS,): row of the first latitude of each N-S band in NSTA/NONL
    11) NSTA (KNREW, KDGL+KNRNS-1): first grid point (1-based) of each row in each E-W
        region; a latitude split between two N-S bands has a row in each
    12) NONL (KNREW, KDGL+KNRNS-1): number of grid points of each row in each E-W region
    13) NGPTOTL (KNREW, KNRNS): grid points of each region
    14) PMU (KDGL,): sines of the latitudes, float64

    Raises RuntimeError with error code 1 when KNRNS or KNREW differs from the region counts
    of the setup, 2 when KDGL differs from the number of latitude rows of KRESOL and 101
    when KRESOL is not a resolution set up by setup_trans_4py.
    """
    return ([KRESOL, KDGL, KNRNS, KNREW],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, (KNRNS,), OUT),
             (np.int64, (KNRNS,), OUT),
             (np.int64, (KNRNS,), OUT),
             (np.int64, (KNREW, KDGL + KNRNS - 1), OUT),
             (np.int64, (KNREW, KDGL + KNRNS - 1), OUT),
             (np.int64, (KNREW, KNRNS), OUT),
             (np.float64, (KDGL,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def trans_inq_vset4py(KRESOL=1):
    """
    Wrapper to TRANS_INQ for the spectral decomposition of a distributed resolution,
    as needed by the *_vset4py bindings.

    Args:\n
    1) KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)

    Returns:\n
    1) KPRTRW: number of wave sets (the number of V-sets is NPROC / KPRTRW)
    2) KMYSETW: this task's wave set (1-based)
    3) KMYSETV: this task's V-set (1-based)
    """
    return ([KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def sp2gp_lam4py(KSIZEI, KSIZEJ,
                   KPHYSICALSIZEI, KPHYSICALSIZEJ,
                   KTRUNCX, KTRUNCY,
                   KNUMMAXRESOL,
                   KSIZE,
                   LGRADIENT,
                   LREORDER,
                   PDELTAX, PDELTAY,
                   PSPEC):
    """
    Transform spectral coefficients into grid-point values.

    Args:\n
    1,2) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    3,4) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    5,6) KTRUNCX, KTRUNCY: troncatures
    7) KNUMMAXRESOL: maximum number of troncatures handled
    8) KSIZE: size of PSPEC
    9) LGRADIENT: gradient computation
    10) LREORDER: reorder spectral coefficients or not
    11,12) PDELTAX,PDELTAY: resolution along x,y axis
    13) PSPEC: spectral coefficient array

    Returns:\n
    1) PGPT: grid-point field
    2) PGPTM: N-S derivative if LGRADIENT
    3) PGPTL: E-W derivative if LGRADIENT
    """
    return ([KSIZEI, KSIZEJ,
             KPHYSICALSIZEI, KPHYSICALSIZEJ,
             KTRUNCX, KTRUNCY,
             KNUMMAXRESOL,
             KSIZE,
             LGRADIENT,
             LREORDER,
             PDELTAX, PDELTAY,
             PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (bool, None, IN),
             (_REAL, None, IN),
             (_REAL, None, IN),
             (_REAL, (KSIZE,), IN),
             (_REAL, (KSIZEI * KSIZEJ,), OUT),
             (_REAL, (KSIZEI * KSIZEJ,), OUT),
             (_REAL, (KSIZEI * KSIZEJ,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gp2sp_lam4py(KSIZE,
                   KSIZEI, KSIZEJ,
                   KPHYSICALSIZEI, KPHYSICALSIZEJ,
                   KTRUNCX, KTRUNCY,
                   KNUMMAXRESOL,
                   PDELTAX, PDELTAY,
                   LREORDER,
                   PGPT):
    """
    Transform grid point values into spectral coefficients.

    Args:\n
    1) KSIZE: size of spectral field
    2,3) KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
    4,5) KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
    6,7) KTRUNCX, KTRUNCY: troncatures
    8) KNUMMAXRESOL: maximum number of troncatures handled
    9,10) PDELTAX, PDELTAY: resolution along x,y axis
    11) LREORDER: reorder spectral coefficients or not
    12) PGPT: grid-point field

    Returns:\n
    1) PSPEC: spectral coefficient array
    """
    return ([KSIZE,
             KSIZEI, KSIZEJ,
             KPHYSICALSIZEI, KPHYSICALSIZEJ,
             KTRUNCX, KTRUNCY,
             KNUMMAXRESOL,
             PDELTAX, PDELTAY,
             LREORDER,
             PGPT],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, None, IN),
             (_REAL, None, IN),
             (bool, None, IN),
             (_REAL, (KSIZEI * KSIZEJ,), IN),
             (_REAL, (KSIZE,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def sp2gp_gauss4py(KSIZEJ,
                     KTRUNC,
                     KNUMMAXRESOL,
                     KGPTOT,
                     KSLOEN,
                     KLOEN,
                     KSIZE,
                     LGRADIENT,
                     LREORDER,
                     PSPEC):
    """
    Transform spectral coefficients into grid-point values.

    Args:\n
    1) KSIZEJ: Number of latitudes
    2) KTRUNC: troncature
    3) KNUMMAXRESOL: maximum number of troncatures handled
    4) KGPTOT: number of grid-points
    5) KSLOEN: Size of KLOEN
    6) KLOEN:
    7) KSIZE: Size of PSPEC
    8) LGRADIENT: compute derivatives
    9) LREORDER: reorder spectral coefficients or not
    10) PSPEC: spectral coefficient array

    Returns:\n
    1) PGPT: grid-point field
    2) PGPTM: N-S derivative if LGRADIENT
    3) PGPTL: E-W derivative if LGRADIENT
    """
    return ([KSIZEJ,
             KTRUNC,
             KNUMMAXRESOL,
             KGPTOT,
             KSLOEN,
             KLOEN,
             KSIZE,
             LGRADIENT,
             LREORDER,
             PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (bool, None, IN),
             (_REAL, (KSIZE,), IN),
             (_REAL, (KGPTOT,), OUT),
             (_REAL, (KGPTOT,), OUT),
             (_REAL, (KGPTOT,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gp2sp_gauss4py(KSPEC,
                     KSIZEJ,
                     KTRUNC,
                     KNUMMAXRESOL,
                     KSLOEN,
                     KLOEN,
                     KSIZE,
                     LREORDER,
                     PGPT):
    """
    Transform grid-point values into spectral coefficients.

    Args:\n
    1) KSPEC: size of spectral coefficients array
    2) KSIZEJ: Number of latitudes
    3) KTRUNC: troncature
    4) KNUMMAXRESOL: maximum number of troncatures handled
    5) KSLOEN: Size of KLOEN
    6) KLOEN
    7) KSIZE: Size of PGPT
    8) LREORDER: reorder spectral coefficients or not
    9) PGPT: grid-point field

    Returns:\n
    1) PSPEC: spectral coefficient array
    """
    return ([KSPEC,
             KSIZEJ,
             KTRUNC,
             KNUMMAXRESOL,
             KSLOEN,
             KLOEN,
             KSIZE,
             LREORDER,
             PGPT],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (_REAL, (KSIZE,), IN),
             (_REAL, (KSPEC,), OUT)],
            None)


@ctypesFF()
def sp2gp_fft1d4py(KSIZES, KTRUNC, PSPEC, KSIZEG):
    """
    Transform spectral coefficients into grid-point values,
    for a 1D array (vertical section academic model)

    Args:\n
    1) KSIZES size of PSPEC
    2) KTRUNC: troncature
    3) PSPEC: spectral coefficient array
    4) KSIZEG: size of grid-point field (with extension zone)

    Returns:\n
    1) PGPT: grid-point field
    """
    return ([KSIZES, KTRUNC, PSPEC, KSIZEG],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KSIZES,), IN),
             (np.int64, None, IN),
             (_REAL, (KSIZEG,), OUT)],
            None)

# ---------------------------------------------------------------------------
# Adjoint distributed spectral transforms on LOCAL arrays
# adj(INV_TRANS) = S^T != S^T W = DIR_TRANS (differ by Gaussian weights W).
# Inner-product convention: <INV_TRANS x, y>_Euclid = <x, INV_TRANSAD y>_mfold
#   where _Euclid is the unweighted grid sum and _mfold applies factor 2 for m>0.
# ---------------------------------------------------------------------------


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_dist4py_ad(KSPEC2, KGPTOT, KFLD, PGP, KRESOL=1):
    """Adjoint of inv_trans_scalar_dist4py (INV_TRANSAD): grid-point seed -> spectral result.
    PGP (KFLD, KGPTOT) IN; returns PSPEC (KFLD, KSPEC2).
    Note: adj(INV_TRANS) = S^T != DIR_TRANS = S^T W (differ by Gaussian weights W).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return (
        [KSPEC2, KGPTOT, KFLD, PGP, KRESOL],
        [
            (np.int64, None, IN),
            (np.int64, None, IN),
            (np.int64, None, IN),
            (_REAL, (KFLD, KGPTOT), IN),
            (_REAL, (KFLD, KSPEC2), OUT),
            (np.int64, None, IN),
        ],
        None,
    )


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_ders_dist4py_ad(KSPEC2, KGPTOT, KFLD, PGP, PGPNS, PGPEW, KRESOL=1):
    """Adjoint of inv_trans_scalar_ders_dist4py (INV_TRANSAD, LDSCDERS=.TRUE.):
    grid value seed PGP, N-S derivative seed PGPNS, E-W derivative seed PGPEW
    (each (KFLD, KGPTOT)) IN; returns the accumulated spectral result PSPEC
    (KFLD, KSPEC2). Exposes the real INV_TRANSAD's own LDSCDERS support (the
    derivative adjoint is already implemented in ecTrans -- this wrapper adds no
    new transform maths, only passes the seeds through in the field order
    inv_trans_scalar_ders_dist4py itself uses: values | N-S | E-W).
    Note: adj(INV_TRANS) = S^T != DIR_TRANS = S^T W (differ by Gaussian weights W).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return (
        [KSPEC2, KGPTOT, KFLD, PGP, PGPNS, PGPEW, KRESOL],
        [
            (np.int64, None, IN),
            (np.int64, None, IN),
            (np.int64, None, IN),
            (_REAL, (KFLD, KGPTOT), IN),
            (_REAL, (KFLD, KGPTOT), IN),
            (_REAL, (KFLD, KGPTOT), IN),
            (_REAL, (KFLD, KSPEC2), OUT),
            (np.int64, None, IN),
        ],
        None,
    )


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_scalar_dist4py_ad(KSPEC2, KGPTOT, KFLD, PSPEC, KRESOL=1):
    """Adjoint of dir_trans_scalar_dist4py (DIR_TRANSAD): spectral seed -> grid result.
    PSPEC (KFLD, KSPEC2) IN; returns PGP (KFLD, KGPTOT).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return (
        [KSPEC2, KGPTOT, KFLD, PSPEC, KRESOL],
        [
            (np.int64, None, IN),
            (np.int64, None, IN),
            (np.int64, None, IN),
            (_REAL, (KFLD, KSPEC2), IN),
            (_REAL, (KFLD, KGPTOT), OUT),
            (np.int64, None, IN),
        ],
        None,
    )


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_uv_dist4py_ad(KSPEC2, KGPTOT, KFLD, PGPU, PGPV, KRESOL=1):
    """Adjoint of inv_trans_uv_dist4py (INV_TRANSAD): grid u,v seeds -> spectral vor,div.
    PGPU, PGPV (KFLD, KGPTOT) IN; returns PSPVOR, PSPDIV (KFLD, KSPEC2).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return (
        [KSPEC2, KGPTOT, KFLD, PGPU, PGPV, KRESOL],
        [
            (np.int64, None, IN),
            (np.int64, None, IN),
            (np.int64, None, IN),
            (_REAL, (KFLD, KGPTOT), IN),
            (_REAL, (KFLD, KGPTOT), IN),
            (_REAL, (KFLD, KSPEC2), OUT),
            (_REAL, (KFLD, KSPEC2), OUT),
            (np.int64, None, IN),
        ],
        None,
    )


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_uv_dist4py_ad(KSPEC2, KGPTOT, KFLD, PSPVOR, PSPDIV, KRESOL=1):
    """Adjoint of dir_trans_uv_dist4py (DIR_TRANSAD): spectral vor,div seeds -> grid u,v.
    PSPVOR, PSPDIV (KFLD, KSPEC2) IN; returns PGPU, PGPV (KFLD, KGPTOT).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return (
        [KSPEC2, KGPTOT, KFLD, PSPVOR, PSPDIV, KRESOL],
        [
            (np.int64, None, IN),
            (np.int64, None, IN),
            (np.int64, None, IN),
            (_REAL, (KFLD, KSPEC2), IN),
            (_REAL, (KFLD, KSPEC2), IN),
            (_REAL, (KFLD, KGPTOT), OUT),
            (_REAL, (KFLD, KGPTOT), OUT),
            (np.int64, None, IN),
        ],
        None,
    )


# === distributed-memory (MPI) interface ===


@ctypesFF()
def mpl_init4py():
    """Initialise FIAT MPL (after mpi4py). Returns (rank[1-based], size)."""
    return ([],
            [(np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@ctypesFF()
def mpl_init_comm4py(KCOMM):
    """Initialise FIAT MPL (after mpi4py) on the communicator KCOMM instead of
    MPI_COMM_WORLD: ecTrans then runs on the tasks of KCOMM only, and the other
    tasks take no part. Call it instead of mpl_init4py, on the tasks of KCOMM.
    KCOMM is the Fortran handle of the communicator (mpi4py: comm.py2f()).
    Returns (rank[1-based], size) within KCOMM."""
    return ([KCOMM],
            [(np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@ctypesFF()
def mpl_end4py():
    """Finalise FIAT MPL (does not finalise MPI)."""
    return ([], [], None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def setup_trans0_4py(KPRGPNS, KPRGPEW, KPRTRW, LDEQ_REGIONS, KMAX_RESOL, LDMPOFF):
    """Parallel resolution-independent setup (processor grid).
    Returns (k_regions_ns, k_regions_ew)."""
    return ([KPRGPNS, KPRGPEW, KPRTRW, LDEQ_REGIONS, KMAX_RESOL, LDMPOFF],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def setup_trans_4py(KSMAX, KDGL, KSLOEN, KLOEN, LDSPLIT, LDUSEFLT, LDGRIDONLY=False):
    """Parallel resolution-dependent setup (LDUSEFLT selects Fast Legendre Transform).
    LDGRIDONLY=True sets up grid-point space only, as the IFS does for its radiation grid
    (KSMAX is then not used; the IFS passes 0): dist_grid4py, gath_grid4py and
    trans_inq_grid4py work on such a resolution, every binding that needs spectral or
    Fourier space (the transforms, the spectral scatter/gather, the norms -- gpnorm_trans4py
    too -- and trans_inq4py, trans_inq_vset4py) raises RuntimeError with error code 102.
    CPU only: the GPU version of ecTrans's SETUP_TRANS resets LDGRIDONLY to False and sets
    up truncation 1 instead (src/trans/gpu/external/setup_trans.F90).
    Returns kresol."""
    return ([KSMAX, KDGL, KSLOEN, KLOEN, LDSPLIT, LDUSEFLT, LDGRIDONLY],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (bool, None, IN),
             (bool, None, IN),
             (np.int64, None, OUT),
             (bool, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def setup_trans0_regions4py(KPRGPNS, KPRGPEW, KPRTRW, LDEQ_REGIONS, KMAX_RESOL, LDMPOFF):
    """setup_trans0_4py that also returns the number of E-W regions of each N-S band of the
    grid-point partition (SETUP_TRANS0's K_REGIONS; with equal regions, it differs between
    the bands). The same setup otherwise: call one of the two, once per process.
    Returns (k_regions_ns, k_regions_ew, k_regions), k_regions of length KPRGPNS * KPRGPEW:
    the number of E-W regions of N-S band 1..k_regions_ns, zero beyond."""
    return ([KPRGPNS, KPRGPEW, KPRTRW, LDEQ_REGIONS, KMAX_RESOL, LDMPOFF],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, (KPRGPNS * KPRGPEW,), OUT)],
            None)


# ---------------------------------------------------------------------------
# Distributed data movement (global on a 'from'/'to' MPL rank <-> local).
# Array convention (ctypesForFortran indexing='C'): Python arrays are field-first,
# i.e. the reverse of the Fortran (X, KFLD) dims -> Python (KFLD, X).
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def dist_spec4py(KSPEC2G, KSPEC2, KFLD, KFROM, PSPECG, KRESOL=1):
    """Scatter global spectral PSPECG(KFLD,KSPEC2G) -> local PSPEC(KFLD,KSPEC2).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2G, KSPEC2, KFLD, KFROM, PSPECG, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KSPEC2G), IN),
             (_REAL, (KFLD, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gath_spec4py(KSPEC2G, KSPEC2, KFLD, KTO, PSPEC, KRESOL=1):
    """Gather local spectral PSPEC(KFLD,KSPEC2) -> global PSPECG(KFLD,KSPEC2G).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2G, KSPEC2, KFLD, KTO, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KSPEC2G), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dist_grid4py(KGPTOTG, KGPTOT, KFLD, KFROM, PGPG, KRESOL=1):
    """Scatter global grid PGPG(KFLD,KGPTOTG) -> local PGP(KFLD,KGPTOT).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KGPTOTG, KGPTOT, KFLD, KFROM, PGPG, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KGPTOTG), IN),
             (_REAL, (KFLD, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gath_grid4py(KGPTOTG, KGPTOT, KFLD, KTO, PGP, KRESOL=1):
    """Gather local grid PGP(KFLD,KGPTOT) -> global PGPG(KFLD,KGPTOTG).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KGPTOTG, KGPTOT, KFLD, KTO, PGP, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KGPTOTG), OUT),
             (np.int64, None, IN)],
            None)


# ---------------------------------------------------------------------------
# Distributed spectral transforms on LOCAL arrays (LDSPLIT resolution). Single
# precision; spectral (KFLD,KSPEC2) model order, grid (KFLD,KGPTOT) single block.
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_dist4py(KSPEC2, KGPTOT, KFLD, PSPEC, KRESOL=1):
    """Local inverse transform: scalar spectral (KFLD,KSPEC2) -> grid (KFLD,KGPTOT).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLD, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_ders_dist4py(KSPEC2, KGPTOT, KFLD, PSPEC, KRESOL=1):
    """Local inverse transform with derivatives: scalar spectral (KFLD,KSPEC2) ->
    grid value, N-S derivative, E-W derivative (each (KFLD,KGPTOT)).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLD, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KGPTOT), OUT),
             (_REAL, (KFLD, KGPTOT), OUT),
             (_REAL, (KFLD, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_scalar_dist4py(KSPEC2, KGPTOT, KFLD, PGP, KRESOL=1):
    """Local direct transform: scalar grid (KFLD,KGPTOT) -> spectral (KFLD,KSPEC2).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLD, PGP, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_uv_dist4py(KSPEC2, KGPTOT, KFLD, PSPVOR, PSPDIV, KRESOL=1):
    """Local inverse transform: vorticity/divergence -> u,v (KFLD,KGPTOT each).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLD, PSPVOR, PSPDIV, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KGPTOT), OUT),
             (_REAL, (KFLD, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_uv_dist4py(KSPEC2, KGPTOT, KFLD, PGPU, PGPV, KRESOL=1):
    """Local direct transform: u,v -> vorticity/divergence (KFLD,KSPEC2 each).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLD, PGPU, PGPV, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KSPEC2), OUT),
             (_REAL, (KFLD, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


# ---------------------------------------------------------------------------
# Norms (global, gathered). Spectral input (KFLD, KSPEC2); grid input (KFLD, KGPTOT).
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def specnorm4py(KSPEC2, KFLD, PSPEC, KRESOL=1):
    """Global spectral L2 norm per field. Returns PNORM(KFLD).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KFLD, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD,), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gpnorm_trans4py(KGPTOT, KFLD, PGP, KRESOL=1):
    """Global grid-point average/min/max per field. Returns (PAVE, PMIN, PMAX).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KGPTOT, KFLD, PGP, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD,), OUT),
             (_REAL, (KFLD,), OUT),
             (_REAL, (KFLD,), OUT),
             (np.int64, None, IN)],
            None)


# ---------------------------------------------------------------------------
# V-set-aware distributed bindings, for the 2-D spectral decomposition (several
# V-sets, KPRTRW < NPROC in setup_trans0_4py): the spectral fields are distributed
# over the V-sets as well. KVSET/KVSETSC/KVSETUV (KFLDG,) give the V-set (1-based)
# of each of the KFLDG global fields; a task holds the spectral coefficients of the
# KFLDL fields of its own V-set (KMYSETV of trans_inq_vset4py), in the order of the
# global fields, and the grid-point values of all KFLDG fields at its grid points.
# With a single V-set (KVSET all 1, KFLDL == KFLDG) they compute the same as the
# *_dist4py bindings, which ecTrans does not accept with several V-sets.
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def dist_spec_vset4py(KSPEC2G, KSPEC2, KFLDG, KFLDL, KFROM, KVSET, PSPECG, KRESOL=1):
    """Scatter global spectral PSPECG(KFLDG,KSPEC2G) from the tasks KFROM (as for
    dist_spec4py) -> local PSPEC(KFLDL,KSPEC2), the fields of this task's V-set.
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2G, KSPEC2, KFLDG, KFLDL, KFROM, KVSET, PSPECG, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDG, KSPEC2G), IN),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gath_spec_vset4py(KSPEC2G, KSPEC2, KFLDG, KFLDL, KTO, KVSET, PSPEC, KRESOL=1):
    """Gather local spectral PSPEC(KFLDL,KSPEC2), the fields of this task's V-set ->
    global PSPECG(KFLDG,KSPEC2G) on the tasks KTO (as for gath_spec4py).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2G, KSPEC2, KFLDG, KFLDL, KTO, KVSET, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDG, KSPEC2G), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_vset4py(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PSPEC, KRESOL=1):
    """Local inverse transform: scalar spectral (KFLDL,KSPEC2) -> grid (KFLDG,KGPTOT).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_ders_vset4py(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PSPEC, KRESOL=1):
    """Local inverse transform with derivatives: scalar spectral (KFLDL,KSPEC2) ->
    grid value, N-S derivative, E-W derivative (each (KFLDG,KGPTOT)).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_scalar_vset4py(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PGP, KRESOL=1):
    """Local direct transform: scalar grid (KFLDG,KGPTOT) -> spectral (KFLDL,KSPEC2).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PGP, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_uv_vset4py(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PSPVOR, PSPDIV, KRESOL=1):
    """Local inverse transform: vorticity/divergence (KFLDL,KSPEC2 each) -> u,v
    (KFLDG,KGPTOT each).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PSPVOR, PSPDIV, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_uv_vset4py(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PGPU, PGPV, KRESOL=1):
    """Local direct transform: u,v (KFLDG,KGPTOT each) -> vorticity/divergence
    (KFLDL,KSPEC2 each).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PGPU, PGPV, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_vset4py_ad(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PGP, KRESOL=1):
    """Adjoint of inv_trans_scalar_vset4py (INV_TRANSAD): grid-point seed PGP
    (KFLDG,KGPTOT) IN; returns PSPEC (KFLDL,KSPEC2).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PGP, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_ders_vset4py_ad(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PGP, PGPNS, PGPEW,
                                     KRESOL=1):
    """Adjoint of inv_trans_scalar_ders_vset4py (INV_TRANSAD, LDSCDERS=.TRUE.): value,
    N-S and E-W derivative seeds (each (KFLDG,KGPTOT)) IN; returns the accumulated
    spectral result PSPEC (KFLDL,KSPEC2).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PGP, PGPNS, PGPEW, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_scalar_vset4py_ad(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PSPEC, KRESOL=1):
    """Adjoint of dir_trans_scalar_vset4py (DIR_TRANSAD): spectral seed PSPEC
    (KFLDL,KSPEC2) IN; returns PGP (KFLDG,KGPTOT).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETSC, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_uv_vset4py_ad(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PGPU, PGPV, KRESOL=1):
    """Adjoint of inv_trans_uv_vset4py (INV_TRANSAD): u,v seeds (KFLDG,KGPTOT each) IN;
    returns PSPVOR, PSPDIV (KFLDL,KSPEC2 each).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PGPU, PGPV, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDG, KGPTOT), IN),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (_REAL, (KFLDL, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_uv_vset4py_ad(KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PSPVOR, PSPDIV, KRESOL=1):
    """Adjoint of dir_trans_uv_vset4py (DIR_TRANSAD): vorticity/divergence seeds
    (KFLDL,KSPEC2 each) IN; returns PGPU, PGPV (KFLDG,KGPTOT each).
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFLDL, KFLDG, KVSETUV, PSPVOR, PSPDIV, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (_REAL, (KFLDG, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def specnorm_vset4py(KSPEC2, KFLDG, KFLDL, KVSET, PSPEC, KRESOL=1):
    """Global spectral L2 norm per field from the local PSPEC(KFLDL,KSPEC2), the fields
    of this task's V-set. Returns PNORM(KFLDG), the norms of all global fields, on
    the first task.
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KFLDG, KFLDL, KVSET, PSPEC, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLDG,), IN),
             (_REAL, (KFLDL, KSPEC2), IN),
             (_REAL, (KFLDG,), OUT),
             (np.int64, None, IN)],
            None)


# ---------------------------------------------------------------------------
# One transform for a whole set of fields: winds and scalars together in a
# single INV_TRANS or DIR_TRANS, as the IFS transforms its model state. The
# V-set maps and field counts are those of the *_vset4py bindings (with a single
# V-set: KVSETUV, KVSETSC all 1, KFUV == KFUVG, KFSC == KFSCG).
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_all_dist4py(KSPEC2, KGPTOT, KFUV, KFSC, KFUVG, KFSCG, KFGP, KDDIVGP, KDSCDERS,
                          KVSETUV, KVSETSC, PSPVOR, PSPDIV, PSPSCALAR, KRESOL=1):
    """Local inverse transform of vorticity/divergence (KFUV,KSPEC2 each) and scalar fields
    (KFSC,KSPEC2), the fields of this task's V-set, in a single INV_TRANS. KVSETUV (KFUVG,)
    and KVSETSC (KFSCG,) give the V-set of each global wind and scalar field. KDDIVGP and
    KDSCDERS (0 or 1) ask for the grid-point divergence and the derivatives of the scalars.
    Returns PGP (KFGP,KGPTOT), all fields in the order of inv_trans.h:
      [divergence (KFUVG, if KDDIVGP)] u (KFUVG) v (KFUVG) scalars (KFSCG)
      [N-S derivatives (KFSCG)] [E-W derivatives (KFSCG)] (the last two if KDSCDERS),
    so KFGP = (3 if KDDIVGP else 2) * KFUVG + (3 if KDSCDERS else 1) * KFSCG.
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFUV, KFSC, KFUVG, KFSCG, KFGP, KDDIVGP, KDSCDERS, KVSETUV,
             KVSETSC, PSPVOR, PSPDIV, PSPSCALAR, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFUVG,), IN),
             (np.int64, (KFSCG,), IN),
             (_REAL, (KFUV, KSPEC2), IN),
             (_REAL, (KFUV, KSPEC2), IN),
             (_REAL, (KFSC, KSPEC2), IN),
             (_REAL, (KFGP, KGPTOT), OUT),
             (np.int64, None, IN)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_all_dist4py(KSPEC2, KGPTOT, KFUV, KFSC, KFUVG, KFSCG, KFGP, KVSETUV, KVSETSC, PGP,
                          KRESOL=1):
    """Local direct transform of grid-point fields PGP (KFGP,KGPTOT), all fields in the order
    of dir_trans.h, u (KFUVG) v (KFUVG) scalars (KFSCG), so KFGP = 2 * KFUVG + KFSCG, in a
    single DIR_TRANS. KVSETUV (KFUVG,) and KVSETSC (KFSCG,) give the V-set of each global wind
    and scalar field. Returns PSPVOR, PSPDIV (KFUV,KSPEC2 each) and PSPSCALAR (KFSC,KSPEC2),
    the fields of this task's V-set.
    KRESOL: resolution tag from setup_trans_4py (default 1, the first resolution)."""
    return ([KSPEC2, KGPTOT, KFUV, KFSC, KFUVG, KFSCG, KFGP, KVSETUV, KVSETSC, PGP, KRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFUVG,), IN),
             (np.int64, (KFSCG,), IN),
             (_REAL, (KFGP, KGPTOT), IN),
             (_REAL, (KFUV, KSPEC2), OUT),
             (_REAL, (KFUV, KSPEC2), OUT),
             (_REAL, (KFSC, KSPEC2), OUT),
             (np.int64, None, IN)],
            None)


# ---------------------------------------------------------------------------
# The repacking of spectral arrays between the layout of the Python callers
# (KFLD, KSPEC2) and that of ecTrans, which every distributed binding with
# spectral fields does (SPEC_LAYOUT4PY_MOD in spec_layout4py.F90): by default
# with array-section loops; with the environment variable ECTRANS4PY_OMP set to a
# value that does not start with '0', with loops blocked over the spectral
# coefficients and threaded with OpenMP (OMP_NUM_THREADS), which pay off with
# several threads per task. The results are bit-identical. ECTRANS4PY_OMP is read
# once per process. These bindings expose the repacking, to test and time it.
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def spec_to_fld4py(KSPEC2, KFLD, PSPEC):
    """The transpose (KSPEC2,KFLD) of PSPEC (KFLD,KSPEC2), as the bindings repack their
    spectral inputs for ecTrans."""
    return ([KSPEC2, KFLD, PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KSPEC2, KFLD), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def fld_to_spec4py(KSPEC2, KFLD, PFLD):
    """The transpose (KFLD,KSPEC2) of PFLD (KSPEC2,KFLD), as the bindings repack the spectral
    outputs of ecTrans."""
    return ([KSPEC2, KFLD, PFLD],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KSPEC2, KFLD), IN),
             (_REAL, (KFLD, KSPEC2), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def spec_layout_inq4py():
    """Inquire the repacking of spectral arrays. Returns:
    1) KBLK: number of spectral coefficients in a block of the blocked loops
    2) KBLOCKED: 1 if ECTRANS4PY_OMP selects the blocked, threaded loops, else 0
    3) KTHREADS: number of OpenMP threads of the blocked loops (1 without OpenMP)
    4) KOPENMP: 1 if the repacking is compiled with OpenMP, else 0
    5) KOPENMP_BUILD: 1 if the build asked for OpenMP (ecTrans's feature OMP), else 0;
       KOPENMP must be the same, or the blocked loops run on a single thread."""
    return ([],
            [(np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


__version__ = ectrans_version().strip()
