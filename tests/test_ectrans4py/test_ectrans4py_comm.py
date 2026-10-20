# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# ectrans4py on a communicator split from MPI_COMM_WORLD: MPL, and with it ecTrans, is
# initialised by mpl_init_comm4py on the communicator of the even ranks, and the odd ranks take
# no part in it. A process initialises MPL once, so this module runs in an MPI job of its own:
#     mpirun -np 4 python -m pytest -k test_ectrans4py_comm

import platform
from types import SimpleNamespace
from unittest import TestCase

import numpy as np
import pytest

import ectrans4py

mpi4py = pytest.importorskip("mpi4py")
from mpi4py import MPI

TRUNCATION = 21
NUM_LATS = 32
LONS_PER_LAT = np.full(NUM_LATS, 2 * NUM_LATS, dtype=np.int64)
NFLD = 4

# Spectral -> grid -> spectral on a regular quadratic Gaussian grid: rounding only
ROUND_TRIP_TOL = 1000 * np.finfo(ectrans4py._REAL).eps

# The communicator ecTrans runs on, and its transform sizes, set up once by setUpModule
context = SimpleNamespace()


def _holds_on_comm(check):
    """Whether check() holds on all tasks of the communicator, returned on every task of
    MPI_COMM_WORLD: the tasks outside the communicator do not call check, and take part only
    in the reduction of the result (with mpi4py)."""
    holds = bool(check()) if context.in_comm else True
    return MPI.COMM_WORLD.allreduce(holds, op=MPI.LAND)


def _root(nfld):
    """MPL rank (1-based) of the task holding the global fields: the first of the communicator"""
    return np.ones(nfld, dtype=np.int64)


def _random_spectral_field(seed):
    """Random global spectral field of real fields: no imaginary part for the zonal wavenumber
    m = 0, whose coefficients come first"""
    rng = np.random.default_rng(seed)
    field = rng.uniform(-1.0, 1.0, (NFLD, context.num_spec2_glob)).astype(ectrans4py._REAL)
    field[:, 1:2 * (TRUNCATION + 1):2] = 0.0
    return field


def setUpModule():
    world = MPI.COMM_WORLD
    context.in_comm = world.Get_rank() % 2 == 0
    context.comm = world.Split(color=0 if context.in_comm else 1, key=world.Get_rank())
    if not context.in_comm:
        return

    system = platform.system()
    if system == "Linux":
        ectrans4py.init_env(unlimited_stack=True)
    elif system == "Darwin":
        ectrans4py.init_env(unlimited_stack=False)
    else:
        raise NotImplementedError("ectrans4py does not support Windows")

    context.mpl_rank, context.mpl_size = ectrans4py.mpl_init4py(context.comm.py2f())
    n_rank = context.comm.Get_size()

    # KPRGPNS=n_rank, KPRGPEW=1, KPRTRW=n_rank, LDEQ_REGIONS=False, KMAX_RESOL=1, LDMPOFF=False
    ectrans4py.setup_trans0_4py(n_rank, 1, n_rank, False, 1, False)

    # KSMAX=TRUNCATION, KDGL=NUM_LATS, KSLOEN=NUM_LATS, KLOEN=LONS_PER_LAT, LDSPLIT=True,
    # LDUSEFLT=False
    kresol = ectrans4py.setup_trans_4py(TRUNCATION, NUM_LATS, NUM_LATS, LONS_PER_LAT, True, False)
    (context.num_grid_points_loc, _, context.num_spec2_loc, context.num_grid_points_glob,
     context.num_spec2_glob, *_) = \
        ectrans4py.trans_inq4py(kresol, NUM_LATS, 0, 0, np.array([], dtype=np.int64), 0)


def tearDownModule():
    if context.in_comm:
        ectrans4py.mpl_end4py()
    context.comm.Free()


class TestCommunicator(TestCase):
    def test_mpl_tasks(self):
        # MPL has the tasks of the communicator only, numbered as in the communicator
        def check():
            return (context.mpl_size == context.comm.Get_size()
                    and context.mpl_rank == context.comm.Get_rank() + 1)

        assert _holds_on_comm(check)

    def test_partition(self):
        # The grid points and spectral coefficients are distributed over the tasks of the
        # communicator only
        def check():
            comm = context.comm
            num_grid_points = comm.allreduce(context.num_grid_points_loc, op=MPI.SUM)
            num_spec2 = comm.allreduce(context.num_spec2_loc, op=MPI.SUM)
            return (context.num_grid_points_glob == np.sum(LONS_PER_LAT) == num_grid_points
                    and context.num_spec2_glob == (TRUNCATION + 1) * (TRUNCATION + 2)
                    == num_spec2)

        assert _holds_on_comm(check)

    def test_round_trip(self):
        # spectral -> grid -> spectral, the global fields on the first task of the communicator
        def check():
            spectral_field_glob = _random_spectral_field(seed=1)
            spectral_field_loc = ectrans4py.dist_spec4py(
                context.num_spec2_glob, context.num_spec2_loc, NFLD, _root(NFLD),
                spectral_field_glob)
            gridded_field_loc = ectrans4py.inv_trans_scalar_dist4py(
                context.num_spec2_loc, context.num_grid_points_loc, NFLD, spectral_field_loc)
            round_trip_loc = ectrans4py.dir_trans_scalar_dist4py(
                context.num_spec2_loc, context.num_grid_points_loc, NFLD, gridded_field_loc)
            round_trip_glob = ectrans4py.gath_spec4py(
                context.num_spec2_glob, context.num_spec2_loc, NFLD, _root(NFLD), round_trip_loc)
            if context.comm.Get_rank() != 0:
                return True
            error = np.max(np.abs(round_trip_glob.astype(np.float64) - spectral_field_glob))
            return error <= ROUND_TRIP_TOL * np.max(np.abs(spectral_field_glob))

        assert _holds_on_comm(check)

    def test_dist_gath(self):
        # Distributing global fields from the first task of the communicator and gathering them
        # back only moves data: bit for bit the same fields
        def check():
            spectral_field_glob = _random_spectral_field(seed=2)
            spectral_field_loc = ectrans4py.dist_spec4py(
                context.num_spec2_glob, context.num_spec2_loc, NFLD, _root(NFLD),
                spectral_field_glob)
            gathered_spectral_glob = ectrans4py.gath_spec4py(
                context.num_spec2_glob, context.num_spec2_loc, NFLD, _root(NFLD),
                spectral_field_loc)
            gridded_field_glob = np.random.default_rng(3).uniform(
                -1.0, 1.0, (NFLD, context.num_grid_points_glob)).astype(ectrans4py._REAL)
            gridded_field_loc = ectrans4py.dist_grid4py(
                context.num_grid_points_glob, context.num_grid_points_loc, NFLD, _root(NFLD),
                gridded_field_glob)
            gathered_gridded_glob = ectrans4py.gath_grid4py(
                context.num_grid_points_glob, context.num_grid_points_loc, NFLD, _root(NFLD),
                gridded_field_loc)
            if context.comm.Get_rank() != 0:
                return True
            return (np.array_equal(gathered_spectral_glob, spectral_field_glob)
                    and np.array_equal(gathered_gridded_glob, gridded_field_glob))

        assert _holds_on_comm(check)
