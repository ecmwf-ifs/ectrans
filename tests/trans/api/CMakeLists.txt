# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Function for generating a suite of api tests for one external subroutine using CMake's
# create_test_sourcelist function
function( generate_api_test_suite )
  set( options )
  set( oneValueArgs SUITE_NAME)
  set( multiValueArgs TEST_LIST WILL_FAIL_LIST )
  cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  set( suite_name ${_PAR_SUITE_NAME} )
  set( test_list ${_PAR_TEST_LIST} )
  set( will_fail_list ${_PAR_WILL_FAIL_LIST} )

  # Add api tests for this suite, with a separate driver executable for each precision
  foreach( precision ${precisions} )
    # Create test suite driver executable encapsulating all tests for this suite at this precision
    create_test_sourcelist( _ ${suite_name}_${precision}.c ${test_list} )
    ecbuild_add_executable( TARGET ${suite_name}_${precision}
                            SOURCES
                                ${suite_name}_${precision}.c
                                ${suite_name}.F90
                                ${CMAKE_SOURCE_DIR}/tests/trans/api/util.F90
                            LIBS trans_${precision} )

    # Prevent the Fortran runtime from trying to provide a main program
    target_link_options( ${suite_name}_${precision} PRIVATE ${NO_FORTRAN_MAIN_FLAG} )

    ectrans_target_fortran_module_directory( TARGET ${suite_name}_${precision}
        MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${suite_name}_${precision} )

    # Add a test for each MPI configuration
    foreach( mpi ${ntasks} )
      # Add each test in the test list for this precision and MPI configuration
      foreach( test ${test_list} )
        ecbuild_add_test( TARGET ${test}_${precision}_mpi${mpi}
                          COMMAND ${suite_name}_${precision}
                          ARGS ${test}
                          MPI ${mpi} )

        # Declare that this test should fail
        if( ${test} IN_LIST will_fail_list )
          set_property( TEST ${test}_${precision}_mpi${mpi} PROPERTY WILL_FAIL TRUE )
          set_property( TEST ${test}_${precision}_mpi${mpi} PROPERTY LABELS "will_fail" )
        endif()
      endforeach()
    endforeach()
  endforeach()
endfunction( generate_api_test_suite )

add_subdirectory( setup_trans0 )
