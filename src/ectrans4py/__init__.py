#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ectrans4py:

A Python interface to spectral transforms from ecTrans, using cTypesForFortran for the Fortran/Python binding.

Two interfaces are provided, sharing one geometry inquiry (``trans_inq4py``) and
one set of Legendre assets (``get_legendre_assets``):

* Serial (single task): the LAM/Gaussian wrappers ``sp2gp_*``/``gp2sp_*``,
  ``etrans_inq4py``, ``get_legendre_assets``, and ``trans_inq4py`` with
  ``KRESOL<=0`` (it self-initialises from the grid parameters).
* Distributed-memory (MPI): initialise with ``mpl_init4py`` (attaches to an
  existing communicator, e.g. mpi4py), set up the processor grid/resolution with
  ``setup_trans0_4py``/``setup_trans_4py`` (``LDSPLIT=.TRUE.``), query the local
  geometry with ``trans_inq4py`` passing the returned ``KRESOL``, transform local
  arrays with the ``*_dist4py`` routines, move data with ``dist_spec``/
  ``gath_spec``/``dist_grid``/``gath_grid``, and compute global norms with
  ``specnorm4py``/``gpnorm_trans4py``.
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
    msg = f'libectrans4py_{{dp/sp}}.{platform_ext} was not found in any of potential locations: {str(lpath)}.'
    msg += 'You can specify a different location using env var LD_LIBRARY_PATH'
    raise FileNotFoundError(msg)

# Floating-point type of the field data (spectral / grid-point) crossing the
# Fortran interface (the JPRB arrays): float32 for the single-precision build,
# float64 for the double. Geometry (Gaussian latitudes/weights, Legendre
# polynomials) and resolution deltas remain double regardless.
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
    KNMENG and PGW are specified across the full globe, pole to pole. PRPNM is specified across the Northern hemisphere only.

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
             (_REAL, (KSLOEN,), OUT),
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
      grid parameters are unused.

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
    8) PMU: sines of the Gaussian latitudes (global, length KSIZEJ)
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
             (_REAL, (KSIZEJ,), OUT),
             (_REAL, (KSIZEJ,), OUT)],
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

__version__ = ectrans_version().strip()


# === distributed-memory (MPI) interface ===


@ctypesFF()
def mpl_init4py():
    """Initialise FIAT MPL (after mpi4py). Returns (rank[1-based], size)."""
    return ([],
            [(np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@ctypesFF()
def mpl_end4py():
    """Finalise FIAT MPL (does not finalise MPI)."""
    return ([], [], None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def setup_trans0_4py(KPRGPNS, KPRGPEW, KPRTRW, LDEQ_REGIONS, KMAX_RESOL):
    """Parallel resolution-independent setup (processor grid). Returns (k_regions_ns, k_regions_ew)."""
    return ([KPRGPNS, KPRGPEW, KPRTRW, LDEQ_REGIONS, KMAX_RESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (bool, None, IN),
             (np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def setup_trans_4py(KSMAX, KDGL, KSLOEN, KLOEN, LDSPLIT, LDUSEFLT):
    """Parallel resolution-dependent setup (LDUSEFLT selects Fast Legendre Transform).
    Returns kresol."""
    return ([KSMAX, KDGL, KSLOEN, KLOEN, LDSPLIT, LDUSEFLT],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (bool, None, IN),
             (bool, None, IN),
             (np.int64, None, OUT)],
            None)


# ---------------------------------------------------------------------------
# Distributed data movement (global on a 'from'/'to' MPL rank <-> local).
# Array convention (ctypesForFortran indexing='C'): Python arrays are field-first,
# i.e. the reverse of the Fortran (X, KFLD) dims -> Python (KFLD, X).
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def dist_spec4py(KSPEC2G, KSPEC2, KFLD, KFROM, PSPECG):
    """Scatter global spectral PSPECG(KFLD,KSPEC2G) -> local PSPEC(KFLD,KSPEC2)."""
    return ([KSPEC2G, KSPEC2, KFLD, KFROM, PSPECG],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KSPEC2G), IN),
             (_REAL, (KFLD, KSPEC2), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gath_spec4py(KSPEC2G, KSPEC2, KFLD, KTO, PSPEC):
    """Gather local spectral PSPEC(KFLD,KSPEC2) -> global PSPECG(KFLD,KSPEC2G)."""
    return ([KSPEC2G, KSPEC2, KFLD, KTO, PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KSPEC2G), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dist_grid4py(KGPTOTG, KGPTOT, KFLD, KFROM, PGPG):
    """Scatter global grid PGPG(KFLD,KGPTOTG) -> local PGP(KFLD,KGPTOT)."""
    return ([KGPTOTG, KGPTOT, KFLD, KFROM, PGPG],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KGPTOTG), IN),
             (_REAL, (KFLD, KGPTOT), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gath_grid4py(KGPTOTG, KGPTOT, KFLD, KTO, PGP):
    """Gather local grid PGP(KFLD,KGPTOT) -> global PGPG(KFLD,KGPTOTG)."""
    return ([KGPTOTG, KGPTOT, KFLD, KTO, PGP],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KFLD,), IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KGPTOTG), OUT)],
            None)


# ---------------------------------------------------------------------------
# Distributed spectral transforms on LOCAL arrays (LDSPLIT resolution). Single
# precision; spectral (KFLD,KSPEC2) model order, grid (KFLD,KGPTOT) single block.
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_dist4py(KSPEC2, KGPTOT, KFLD, PSPEC):
    """Local inverse transform: scalar spectral (KFLD,KSPEC2) -> grid (KFLD,KGPTOT)."""
    return ([KSPEC2, KGPTOT, KFLD, PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KGPTOT), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_scalar_ders_dist4py(KSPEC2, KGPTOT, KFLD, PSPEC):
    """Local inverse transform with derivatives: scalar spectral (KFLD,KSPEC2) ->
    grid value, N-S derivative, E-W derivative (each (KFLD,KGPTOT))."""
    return ([KSPEC2, KGPTOT, KFLD, PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KGPTOT), OUT),
             (_REAL, (KFLD, KGPTOT), OUT),
             (_REAL, (KFLD, KGPTOT), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_scalar_dist4py(KSPEC2, KGPTOT, KFLD, PGP):
    """Local direct transform: scalar grid (KFLD,KGPTOT) -> spectral (KFLD,KSPEC2)."""
    return ([KSPEC2, KGPTOT, KFLD, PGP],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KSPEC2), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def inv_trans_uv_dist4py(KSPEC2, KGPTOT, KFLD, PSPVOR, PSPDIV):
    """Local inverse transform: vorticity/divergence -> u,v (KFLD,KGPTOT each)."""
    return ([KSPEC2, KGPTOT, KFLD, PSPVOR, PSPDIV],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD, KGPTOT), OUT),
             (_REAL, (KFLD, KGPTOT), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def dir_trans_uv_dist4py(KSPEC2, KGPTOT, KFLD, PGPU, PGPV):
    """Local direct transform: u,v -> vorticity/divergence (KFLD,KSPEC2 each)."""
    return ([KSPEC2, KGPTOT, KFLD, PGPU, PGPV],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD, KSPEC2), OUT),
             (_REAL, (KFLD, KSPEC2), OUT)],
            None)


# ---------------------------------------------------------------------------
# Norms (global, gathered). Spectral input (KFLD, KSPEC2); grid input (KFLD, KGPTOT).
# ---------------------------------------------------------------------------

@treatReturnCode
@ctypesFF()
@addReturnCode
def specnorm4py(KSPEC2, KFLD, PSPEC):
    """Global spectral L2 norm per field. Returns PNORM(KFLD)."""
    return ([KSPEC2, KFLD, PSPEC],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KSPEC2), IN),
             (_REAL, (KFLD,), OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def gpnorm_trans4py(KGPTOT, KFLD, PGP):
    """Global grid-point average/min/max per field. Returns (PAVE, PMIN, PMAX)."""
    return ([KGPTOT, KFLD, PGP],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (_REAL, (KFLD, KGPTOT), IN),
             (_REAL, (KFLD,), OUT),
             (_REAL, (KFLD,), OUT),
             (_REAL, (KFLD,), OUT)],
            None)


__version__ = ectrans_version().strip()
