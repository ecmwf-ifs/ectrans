# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

from unittest import TestCase
import ectrans4py
import platform
import numpy as np
import pytest
from types import SimpleNamespace

mpi4py = pytest.importorskip("mpi4py")
from mpi4py import MPI

TRUNCATION = 79
NUM_LATS = 2 * (TRUNCATION + 1)
LONS_PER_LAT = [20 + 4 * i for i in range(NUM_LATS // 2)]
LONS_PER_LAT = np.array(LONS_PER_LAT + LONS_PER_LAT[::-1])
NFLD = 20

# Shared transform geometry, populated once by setUpModule
context = SimpleNamespace()

def _from_which_rank(rank, nfld):
    return np.ones(nfld, dtype=np.int64)

def _initialise_global_spectral_field(rank, nfld, num_spec2_glob):
    # Only rank 0 will have non-zero, meaningful data in the global spectral field
    global_spectral_field = np.zeros((nfld, num_spec2_glob), dtype=ectrans4py._REAL)
    if rank == 0:
        global_spectral_field[:, 0] = 1.0 # (m = 0, n = 0) real part is 1.0
    return global_spectral_field

def _initialise_local_spectral_field(rank, nfld, num_spec2_glob, num_spec2_loc):
    global_spectral_field = _initialise_global_spectral_field(rank, nfld, num_spec2_glob)
    return ectrans4py.dist_spec4py(num_spec2_glob, num_spec2_loc, nfld,
                                   _from_which_rank(rank, nfld), global_spectral_field)

def setUpModule():
    system = platform.system()
    if system == "Linux":
        ectrans4py.init_env(unlimited_stack=True)
    elif system == "Darwin":
        ectrans4py.init_env(unlimited_stack=False)
    else:
        raise NotImplementedError("ectrans4py does not support Windows")

    comm = MPI.COMM_WORLD
    n_rank = comm.Get_size()
    my_rank = comm.Get_rank()
    if comm.Get_rank() == 0:
        print(f"Testing with {comm.Get_size()} MPI processes")
    ectrans4py.mpl_init4py()

    # KPRGPNS=n_rank, KPRGPEW=1, KPRTRW=n_rank, LDEQ_REGIONS=False, KMAX_RESOL=1, LDMPOFF=False
    ectrans4py.setup_trans0_4py(n_rank, 1, n_rank, False, 1, False)

    # KSMAX=TRUNCATION, KDGL=num_lats, KSLOEN=num_lats, KLOEN=lons_per_lat, LDSPLIT=True
    # LDUSEFLT=False
    kresol = ectrans4py.setup_trans_4py(TRUNCATION, NUM_LATS, NUM_LATS, LONS_PER_LAT, True, False)

    (num_grid_points_loc, _, num_spec2_loc, num_grid_points_glob, num_spec2_glob, _, _, _, _) = \
        ectrans4py.trans_inq4py(kresol, NUM_LATS, 0, 0, np.array([], dtype=np.int64), 0)

    print(f"num_grid_points_loc: {num_grid_points_loc}, num_spec2_loc: {num_spec2_loc}, num_grid_points_glob: {num_grid_points_glob}, num_spec2_glob: {num_spec2_glob}")

    context.my_rank = my_rank
    context.n_rank = n_rank
    context.kresol = kresol
    context.num_grid_points_loc = num_grid_points_loc
    context.num_spec2_loc = num_spec2_loc
    context.num_grid_points_glob = num_grid_points_glob
    context.num_spec2_glob = num_spec2_glob

def tearDownModule():
    ectrans4py.mpl_end4py()


class TestInverse(TestCase):
    def test_scalar(self):
        local_spectral_field = _initialise_local_spectral_field(context.my_rank, NFLD,
                                                               context.num_spec2_glob,
                                                               context.num_spec2_loc)
        gridded_local = ectrans4py.inv_trans_scalar_dist4py(context.num_spec2_loc,
                                                            context.num_grid_points_loc, NFLD, local_spectral_field)
        gridded_glob = ectrans4py.gath_grid4py(context.num_grid_points_glob,
                                               context.num_grid_points_loc, NFLD,
                                               _from_which_rank(context.my_rank, NFLD),
                                               gridded_local)
        if context.my_rank == 0:
            range_val = float(np.max(gridded_glob) - np.min(gridded_glob))
            print(range_val)
            assert range_val <= 1.0e-4

    def test_scalar_ders(self):
        assert True

    def test_uv(self):
        assert True

class TestDirect(TestCase):
    def test_scalar(self):
        assert True

    def test_uv(self):
        assert True

class TestRoundTrip(TestCase):
    def test_scalar(self):
        assert True

    def test_uv(self):
        assert True

class TestDistributeGather(TestCase):
    def test_dist_spec(self):
        assert True

    def test_dist_grid(self):
        assert True

class TestNorm(TestCase):
    def test_specnorm(self):
        assert True

    def test_gpnorm_trans(self):
        assert True
