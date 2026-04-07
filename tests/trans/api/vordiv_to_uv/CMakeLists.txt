# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# API tests for VORDIV_TO_UV

# Define list of VORDIV_TO_UV tests
set( test_list
  ectrans_test_trans_api_vordiv_to_uv_t1
)

# Declare tests that WILL_FAIL
set( will_fail_list
)

generate_api_test_suite( SUITE_NAME vordiv_to_uv_test_suite
                         TESTS ${test_list}
                         WILL_FAIL_LIST ${will_fail_list}
                         MPIxOMPS 0x1 0x8
                         BACKENDS cpu_sp cpu_dp
                         EXCLUDES  )
# The backends gpu_sp gpu_dp are not yet supported in vordiv_to_uv
