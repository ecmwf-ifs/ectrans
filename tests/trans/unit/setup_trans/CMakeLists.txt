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

# Add unit tests for SETUP_TRANS, with a separate driver executable for each precision
foreach( precision ${precisions} )
  # Create test suite driver executable encapsulating all SETUP_TRANS tests for this precision
  create_test_sourcelist( _ setup_trans_test_suite_${precision}.c ${test_list} )
  ecbuild_add_executable( TARGET setup_trans_test_suite_${precision}
                          SOURCES
                            setup_trans_test_suite_${precision}.c
                            setup_trans_test_suite.F90
                          LIBS trans_${precision} )
  target_link_options( setup_trans_test_suite_${precision} PRIVATE ${NO_FORTRAN_MAIN_FLAG} )
  ectrans_target_fortran_module_directory( TARGET setup_trans_test_suite_${precision}
    MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/setup_trans_test_suite_${precision} )

  # Add each test in the test list for this precision
  foreach( test ${test_list} )
    add_test( NAME ${test}_${precision} COMMAND setup_trans_test_suite_${precision} ${test} )

    # Declare that this test should fail
    if( ${test} IN_LIST will_fail_list )
      set_property( TEST ${test}_${precision} PROPERTY WILL_FAIL TRUE )
      set_property( TEST ${test}_${precision} PROPERTY LABELS "will_fail" )
    endif()
  endforeach()
endforeach()
