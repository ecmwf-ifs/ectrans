# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Check of the repacking of spectral arrays between the layout of the Python callers and that of
ecTrans, which every distributed binding with spectral fields does: a pure permutation, so exact,
with the array-section loops (the default) and with the blocked, OpenMP-threaded loops
(ECTRANS4PY_OMP set), whatever the number of threads. test_ectrans4py_serial.py runs it in its own
process and, as a script that prints its result as JSON, in processes with other settings of
ECTRANS4PY_OMP and OMP_NUM_THREADS, which are read once per process."""

import json

import numpy as np

import ectrans4py

# (number of spectral coefficients, number of fields): below, at and across the size of the blocks
# of the blocked loops (64), several blocks with a remainder, one field, many fields, no field
SHAPES = ((1, 1), (2, 137), (63, 1), (64, 137), (65, 137), (127, 3), (128, 137), (129, 138),
          (1000, 138), (4099, 5), (65, 0))


def report(shapes=SHAPES):
    """The shapes for which the repacking is not exact (both directions, and back), and the
    settings of the repacking (spec_layout_inq4py)"""
    failed = []
    for num_spec2, num_fields in shapes:
        # Values that identify their coefficient and field, exact in single precision
        coefficient = np.arange(num_spec2) % 4096
        field = np.arange(num_fields)[:, np.newaxis] % 256
        spectral = (coefficient + 4096 * field).astype(ectrans4py._REAL)
        fields_first = ectrans4py.spec_to_fld4py(num_spec2, num_fields, spectral)
        back = ectrans4py.fld_to_spec4py(num_spec2, num_fields, fields_first)
        if not (np.array_equal(fields_first, spectral.T) and np.array_equal(back, spectral)):
            failed.append([num_spec2, num_fields])
    block, blocked, threads, openmp, openmp_build = ectrans4py.spec_layout_inq4py()
    return dict(failed=failed, block=int(block), blocked=int(blocked), threads=int(threads),
                openmp=int(openmp), openmp_build=int(openmp_build))


if __name__ == "__main__":
    print(json.dumps(report()))
