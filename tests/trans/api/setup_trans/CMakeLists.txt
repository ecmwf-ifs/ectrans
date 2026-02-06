# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# API tests for SETUP_TRANS

# Define list of SETUP_TRANS tests
set( test_list
  ectrans_test_trans_api_setup_trans_without_setup_trans0
  ectrans_test_trans_api_setup_trans_basic
  ectrans_test_trans_api_setup_trans_multiple_resolutions
  ectrans_test_trans_api_setup_trans_odd_ndgl
  ectrans_test_trans_api_setup_trans_octahedral
  ectrans_test_trans_api_setup_trans_ldsplit
  ectrans_test_trans_api_setup_trans_stretching
  ectrans_test_trans_api_setup_trans_flt
  ectrans_test_trans_api_setup_trans_all_fftw
  ectrans_test_trans_api_setup_trans_belusov
  ectrans_test_trans_api_setup_trans_gridonly
  ectrans_test_trans_api_setup_trans_spsetuponly
  ectrans_test_trans_api_setup_trans_pnmonly
)

# Declare tests that WILL_FAIL
set( will_fail_list
  ectrans_test_trans_api_setup_trans_odd_ndgl
  ectrans_test_trans_api_setup_trans_without_setup_trans0
)

# Note: the FLT is not implemented for the GPU backend
generate_api_test_suite( SUITE_NAME setup_trans_test_suite
                         TESTS ${test_list}
                         WILL_FAIL_LIST ${will_fail_list}
                         MPIxOMPS "0x1" "0x8" "1x1" "1x8" "2x1" "2x8"
                         BACKENDS "cpu_sp" "cpu_dp" "gpu_sp" "gpu_dp"
                         EXCLUDES "flt.*gpu" )
