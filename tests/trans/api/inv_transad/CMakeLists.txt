# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( test_list
  ectrans_test_trans_api_inv_transad_grad_adjoint_t7
  ectrans_test_trans_api_inv_transad_grad_adjoint_t11
  ectrans_test_trans_api_inv_transad_grad_adjoint_t14
)

generate_api_test_suite( SUITE_NAME inv_transad_test_suite
                         TESTS ${test_list}
                         MPIxOMPS 0x1 0x8 1x1 1x8 2x1 2x8
                         BACKENDS cpu_sp cpu_dp gpu_sp gpu_dp )
