# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Function for generating a suite of API tests for one external subroutine using CMake's
# create_test_sourcelist function
function( generate_api_test_suite )
  set( options )
  set( oneValueArgs SUITE_NAME )
  set( multiValueArgs TESTS WILL_FAIL_LIST MPIxOMPS BACKENDS EXCLUDES )
  cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  set( suite_name ${_PAR_SUITE_NAME} )
  set( tests ${_PAR_TESTS} )
  set( will_fail_list ${_PAR_WILL_FAIL_LIST} )
  set( mpixomps ${_PAR_MPIxOMPS} )
  set( backends ${_PAR_BACKENDS} )
  set( excludes ${_PAR_EXCLUDES} )

  # Add API tests for this suite, with a separate driver executable for each backend
  foreach( backend ${backends} )
    # Don't include tests which aren't supported by this build
    if( (backend MATCHES "cpu" AND NOT HAVE_CPU) OR (backend MATCHES "gpu" AND NOT HAVE_GPU) )
      continue()
    endif()
    if( backend MATCHES "sp" AND NOT HAVE_SINGLE_PRECISION )
      continue()
    endif()
    if( backend MATCHES "dp" AND NOT HAVE_DOUBLE_PRECISION )
      continue()
    endif()

    # SETUP_TRANS0 is precision and platform-agnostic
    if( suite_name STREQUAL setup_trans0_test_suite )
      set( libs ectrans_common )
      set( suffix "" )
    else()
      if( backend MATCHES "cpu" )
        set( platform_tag "" ) # The CPU library doesn't have a suffix
      elseif( backend MATCHES "gpu" )
        set( platform_tag "_gpu" )
      endif()
      if( backend MATCHES "sp" )
        set( precision "sp" )
      elseif( backend MATCHES "dp" )
        set( precision "dp" )
      endif()
      set( libs trans${platform_tag}_${precision} )
      set( suffix "_${backend}" )
    endif()

    # Create test suite driver executable encapsulating all tests for this suite at this precision
    create_test_sourcelist( _ ${suite_name}${suffix}.c ${tests} )
    ecbuild_add_executable( TARGET ${suite_name}${suffix}
                            SOURCES
                                ${suite_name}${suffix}.c
                                ${suite_name}.F90
                                ${CMAKE_CURRENT_SOURCE_DIR}/../util.F90
                            LINKER_LANGUAGE Fortran
                            LIBS ${libs} )

    # Prevent the Fortran runtime from trying to provide a main program
    target_link_options( ${suite_name}${suffix} PRIVATE ${NO_FORTRAN_MAIN_FLAG} )

    ecbuild_target_fortran_module_directory( TARGET ${suite_name}${suffix}
        MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${suite_name}${suffix} )

    # Add a test for each MPI x OpenMP configuration
    foreach( mpixomp ${mpixomps} )
      string(REGEX MATCH "^([0-9]+)x([0-9]+)$" _ "${mpixomp}")
      set(mpi "${CMAKE_MATCH_1}")
      set(omp "${CMAKE_MATCH_2}")

      if( (mpi GREATER 0 AND NOT HAVE_MPI) OR (omp GREATER 1 AND NOT HAVE_OMP) )
        continue()
      endif()

      # Add each test in the test list for this precision and MPI configuration
      foreach( test ${tests} )
        # We don't yet have a good way to handle aborts in ecTrans gracefully, so for now we will
        # not include will-fail tests
        if( ${test} IN_LIST will_fail_list )
          continue()
        endif()

        set( full_test_name ${test}${suffix}_mpi${mpi}xomp${omp} )

        # If this test is in the exclusion list, skip it
        set( skip_test FALSE )
        foreach( exclude ${excludes} )
          if( full_test_name MATCHES ${exclude} )
            set( skip_test TRUE )
          endif()
        endforeach()
        if( skip_test )
          continue()
        endif()

        ecbuild_add_test( TARGET ${full_test_name}
                          COMMAND ${suite_name}${suffix}
                          ARGS ${test}
                          MPI ${mpi}
                          OMP ${omp} )

        # Declare that this test should fail
        if( ${test} IN_LIST will_fail_list )
          set_property( TEST ${full_test_name} PROPERTY WILL_FAIL TRUE )
          set_property( TEST ${full_test_name} PROPERTY LABELS "will_fail" )
        endif()
      endforeach()
    endforeach()
  endforeach()
endfunction( generate_api_test_suite )

add_subdirectory( setup_trans0 )
add_subdirectory( setup_trans )
