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

system = platform.system()
if system == "Linux":
    ectrans4py.init_env(unlimited_stack=True)
elif system == "Darwin":
    ectrans4py.init_env(unlimited_stack=False)
else:
    raise NotImplementedError("ectrans4py does not support Windows")

class TestInverse(TestCase):
    def test_scalar(self):
        return True

    def test_scalar_ders(self):
        return True

    def test_uv(self):
        return True

class TestDirect(TestCase):
    def test_scalar(self):
        return True

    def test_uv(self):
        return True

class Distribute(TestCase):
    def dist_spec(self):
        return True

    def dist_grid(self):
        return True

class Gather(TestCase):
    def gath_spec(self):
        return True

    def gath_grid(self):
        return True

class Norm(TestCase):
    def specnorm(self):
        return True

    def gpnorm_trans(self):
        return True
