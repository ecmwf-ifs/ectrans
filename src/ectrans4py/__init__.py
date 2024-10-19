#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
"""
ialsptrans4py:

Contains the interface to spectral transforms from the IAL/ecTrans.
Note that this is temporary between the former package arpifs4py and a direct python interface to ecTrans.

Actual .so library should be in one of the preinstalled paths or in a directory specified via LD_LIBRARY_PATH
"""

from __future__ import print_function, absolute_import, unicode_literals, division

import os
import numpy as np
import ctypesForFortran
from ctypesForFortran import addReturnCode, treatReturnCode, IN, OUT



__version__ = "2.0.1"

# Shared objects library
########################
shared_objects_library = os.environ.get('IALSPTRANS4PY_SO', None)
if shared_objects_library is None or not os.path.exists(shared_objects_library):
    # not specified or path does not exist : find in known locations
    so_basename = "libtrans_dp.so"  # local name in the directory
    LD_LIBRARY_PATH = [p for p in os.environ.get('LD_LIBRARY_PATH', '').split(':') if p != '']
    potential_locations = LD_LIBRARY_PATH + [
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lib'),   # FIXEME : but requiere changes in CMakeLists.txt
        os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lib64'),     # force one libdir directory name ! 
#        "/home/common/epygram/public/EPyGrAM/libs4py",  # CNRM
#        "/home/gmap/mrpe/mary/public/EPyGrAM/libs4py",  # belenos/taranis
#        "/home/acrd/public/EPyGrAM/libs4py",  # ECMWF's Atos aa-ad
        ]
    for _libs4py_dir in potential_locations:
        shared_objects_library = os.path.join(_libs4py_dir, so_basename)
        if os.path.exists(shared_objects_library):
            break
        else:
            shared_objects_library = None
    if shared_objects_library is None:
        msg = ' '.join(["'{}' was not found in any of potential locations: {}.",
                        "You can specify a different location using env var LD_LIBRARY_PATH",
                        "or specify a precise full path using env var IALSPTRANS4PY_SO."]).format(
                                so_filename, str(potential_locations))
        raise FileNotFoundError(msg)
else:
    so_basename = os.path.basename(shared_objects_library)
ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory(shared_objects_library)

# Initialization
################

def init_env(omp_num_threads=None,
             no_mpi=False):
    """
    Set adequate environment for the inner libraries.

    :param int omp_num_threads: sets OMP_NUM_THREADS
    :param bool no_mpi: environment variable DR_HOOK_NOT_MPI set to 1
    """
    # because arpifs library is compiled with MPI & openMP
    if omp_num_threads is not None:
        os.environ['OMP_NUM_THREADS'] = str(omp_num_threads)
    if no_mpi:
        os.environ['DR_HOOK_NOT_MPI'] = '1'

# Transforms interfaces
#######################

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
             (np.float64, None, IN),
             (np.float64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT)],
            None)


@treatReturnCode
@ctypesFF()
@addReturnCode
def trans_inq4py(KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL):
    """
    Simplified wrapper to TRANS_INQ.

    Args:\n
    1) KSIZEJ: number of latitudes in grid-point space
    2) KTRUNC: troncature
    3) KSLOEN: Size of KLOEN
    4) KLOEN: number of points on each latitude row
    5) KNUMMAXRESOL: maximum number of troncatures handled

    Returns:\n
    1) KGPTOT: number of gridpoints
    2) KSPEC: number of spectral coefficients
    3) KNMENG: cut-off zonal wavenumber
    """
    return ([KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL],
            [(np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, None, IN),
             (np.int64, (KSLOEN,), IN),
             (np.int64, None, IN),
             (np.int64, None, OUT),
             (np.int64, None, OUT),
             (np.int64, (KSLOEN,), OUT)],
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
             (np.float64, None, IN),
             (np.float64, None, IN),
             (np.float64, (KSIZE,), IN),
             (np.float64, (KSIZEI * KSIZEJ,), OUT),
             (np.float64, (KSIZEI * KSIZEJ,), OUT),
             (np.float64, (KSIZEI * KSIZEJ,), OUT)],
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
             (np.float64, None, IN),
             (np.float64, None, IN),
             (bool, None, IN),
             (np.float64, (KSIZEI * KSIZEJ,), IN),
             (np.float64, (KSIZE,), OUT)],
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
             (np.float64, (KSIZE,), IN),
             (np.float64, (KGPTOT,), OUT),
             (np.float64, (KGPTOT,), OUT),
             (np.float64, (KGPTOT,), OUT)],
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
             (np.float64, (KSIZE,), IN),
             (np.float64, (KSPEC,), OUT)],
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
             (np.float64, (KSIZES,), IN),
             (np.int64, None, IN),
             (np.float64, (KSIZEG,), OUT)],
            None)
