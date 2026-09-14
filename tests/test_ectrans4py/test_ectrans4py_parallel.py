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
from mpi4py import MPI


class Ectrans4pyTestBase(TestCase):
    _env_initialised = False

    @classmethod
    def setUpClass(cls):
        if not Ectrans4pyTestBase._env_initialised:
            system = platform.system()
            if system == "Linux":
                ectrans4py.init_env(unlimited_stack=True)
            elif system == "Darwin":
                ectrans4py.init_env(unlimited_stack=False)
            else:
                raise NotImplementedError("ectrans4py does not support Windows")
            Ectrans4pyTestBase._env_initialised = True

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        if rank == 0:
            print(f"Testing with {size} MPI processes")
        ectrans4py.mpl_init4py()

    @classmethod
    def tearDownClass(cls):
        ectrans4py.mpl_end4py()


class TestInverse(Ectrans4pyTestBase):
    def test_scalar(self):
        return True

    def test_scalar_ders(self):
        return True

    def test_uv(self):
        return True

class TestDirect(Ectrans4pyTestBase):
    def test_scalar(self):
        return True

    def test_uv(self):
        return True

class Distribute(Ectrans4pyTestBase):
    def dist_spec(self):
        return True

    def dist_grid(self):
        return True

class Gather(Ectrans4pyTestBase):
    def gath_spec(self):
        return True

    def gath_grid(self):
        return True

class Norm(Ectrans4pyTestBase):
    def specnorm(self):
        return True

    def gpnorm_trans(self):
        return True
