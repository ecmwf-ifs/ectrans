# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Unit tests for GPNORM_TRANS

# Define list of GPNORM_TRANS tests
set( test_list
  unit_test_gpnorm_trans_min
  unit_test_gpnorm_trans_max
)

# Declare tests that WILL_FAIL
set( will_fail_list
)

generate_unit_test_suite( SUITE_NAME gpnorm_trans_test_suite
                          TEST_LIST ${test_list}
                          WILL_FAIL_LIST ${will_fail_list} )
