# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Unit tests for SETUP_TRANS

# Define list of SETUP_TRANS tests
set( test_list
  unit_test_setup_trans_without_setup_trans0
  unit_test_setup_trans_basic
  unit_test_setup_trans_odd_ndgl
  unit_test_setup_trans_octahedral
  unit_test_setup_trans_stretching
  unit_test_setup_trans_all_fftw
  unit_test_setup_trans_belusov
  unit_test_setup_trans_flt
)

# Declare tests that WILL_FAIL
set( will_fail_list
  unit_test_setup_trans_odd_ndgl
  unit_test_setup_trans_without_setup_trans0
)

generate_unit_test_suite( SUITE_NAME setup_trans_test_suite
                          TEST_LIST ${test_list}
                          WILL_FAIL_LIST ${will_fail_list} )
