# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# API tests for SETUP_TRANS0

# Define list of SETUP_TRANS0 tests
set( test_list
  ectrans_test_trans_api_setup_trans0_eq_regions
)

# Define additional list setting tests that WILL_FAIL
set( will_fail_list
  "" # None yet
)

generate_api_test_suite( SUITE_NAME setup_trans0_test_suite
                         TESTS ${test_list}
                         WILL_FAIL_LIST ${will_fail_list}
                         MPIxOMPS "0x1" "0x8" "1x1" "1x8" "2x1" "2x8"
                         BACKENDS _ )
