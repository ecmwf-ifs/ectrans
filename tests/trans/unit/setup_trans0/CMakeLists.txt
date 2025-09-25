# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Unit tests for SETUP_TRANS0

# Establish which precisions are enabled
set( precisions "" )
if( HAVE_DOUBLE_PRECISION )
  list( APPEND precisions "dp" )
endif()
if( HAVE_SINGLE_PRECISION )
  list( APPEND precisions "sp" )
endif()

# Define list of SETUP_TRANS0 tests
set( test_list
  unit_test_setup_trans0_eq_regions
)

# Define additional list setting tests that WILL_FAIL
set( will_fail_list
  "" # None yet
)

# Add unit tests for SETUP_TRANS0, with a separate driver executable for each precision
foreach( precision ${precisions} )
  # Create test suite driver executable encapsulating all SETUP_TRANS0 tests for this precision
  create_test_sourcelist( _ setup_trans0_test_suite_${precision}.c ${test_list} )
  ecbuild_add_executable( TARGET setup_trans0_test_suite_${precision}
                          SOURCES
                            setup_trans0_test_suite_${precision}.c
                            setup_trans0_test_suite.F90
                          LIBS trans_${precision} )

  ectrans_target_fortran_module_directory( TARGET setup_trans0_test_suite_${precision}
    MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/setup_trans0_test_suite_${precision} )

  # Add each test in the test list for this precision
  foreach( test ${test_list} )
    add_test( NAME ${test}_${precision} COMMAND setup_trans0_test_suite_${precision} ${test} )

    # Declare that this test should fail
    if( ${test} IN_LIST will_fail_list )
      set_property( TEST ${test}_${precision} PROPERTY WILL_FAIL TRUE )
    endif()
  endforeach()
endforeach()
