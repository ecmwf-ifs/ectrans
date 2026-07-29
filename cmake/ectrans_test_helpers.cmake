# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Test helper functions and macros for ectrans CTest definitions.

# Build a CTest fixture/property-safe name from a logical prefix and target.
# Use this whenever a test name or reference name has to become a fixture or
# global property key, because CTest fixtures are easier to reason about when
# non-identifier characters have been normalized consistently.
function(ectrans_test_fixture_name output prefix target)
  set( fixture_name "${prefix}_${target}" )
  string( REGEX REPLACE "[^A-Za-z0-9_]" "_" fixture_name "${fixture_name}" )
  set( ${output} "${fixture_name}" PARENT_SCOPE )
endfunction()

# Build the global-property base name used to remember a local reference source.
# Use this for canonical checksum references that may be compared by many tests;
# the associated properties record which source test produces the symlink target.
function(ectrans_local_reference_property_name output reference_name)
  ectrans_test_fixture_name( property_name local_reference_file ${reference_name} )
  set( ${output} "${property_name}" PARENT_SCOPE )
endfunction()

# Register a build-local canonical checksum reference.
# Use this for the selected source test that defines the bit-identical baseline
# for a logical checksum case. The function creates a symlink in
# ${local_reference_dir} and records the producing CTest target/fixture so later
# local-reference diff tests can pull the required source test into filtered
# CTest runs. The checksum file prefix basename must match the CTest target that
# generates it.
function(ectrans_add_local_reference_file test_prefix reference_name)
  get_filename_component( source_target "${test_prefix}" NAME )
  set( local_reference_file "${local_reference_dir}/${reference_name}.checksums" )
  ectrans_test_fixture_name( source_fixture checksum_file ${source_target} )
  ectrans_local_reference_property_name( local_reference_property ${reference_name} )
  file( REMOVE "${local_reference_file}" )
  file( CREATE_LINK "${test_prefix}.checksums" "${local_reference_file}" SYMBOLIC )
  set_property( GLOBAL PROPERTY "${local_reference_property}_SOURCE_TARGET" "${source_target}" )
  set_property( GLOBAL PROPERTY "${local_reference_property}_SOURCE_FIXTURE" "${source_fixture}" )
  set_property( GLOBAL APPEND PROPERTY ECTRANS_LOCAL_REFERENCE_SOURCE_TARGETS "${source_target}" )
endfunction()

# Compose the canonical checksum reference basename for a logical benchmark case.
# Use this before adding external or local checksum diff tests so all comparison
# paths and generated local-reference symlinks share the same stable naming
# convention, independent of MPI/OpenMP/callmode variants.
function(ectrans_checksum_reference_name output benchmark truncation grid nfld nlev niter scders uvders vordiv flt)
  set( ${output}
  "${benchmark}_${truncation}_${grid}_${nfld}_${nlev}_${scders}_${uvders}_${vordiv}_${flt}_${niter}"
       PARENT_SCOPE )
endfunction()

# Add checksum comparison tests for one benchmark-producing CTest target.
# Use this immediately after defining a benchmark test that writes
# ${checksums_dir}/${target}.checksums. It adds the external-reference diff test
# when HAVE_REFERENCE_TESTS is enabled, and adds a bit-identical test when
# HAVE_BITID_TESTS is enabled and ectrans_add_local_reference_file()
# registered a source for the same reference name. Fixtures make filtered CTest
# invocations run the checksum-producing prerequisites before the diff tests.
function(ectrans_add_checksum_reference_tests target reference_name)
  if( NOT HAVE_REFERENCE_TESTS AND NOT HAVE_BITID_TESTS )
    return()
  endif()

  set( checksum_file "${checksums_dir}/${target}.checksums" )
  set( reference_file "${reference_name}.checksums" )
  set( local_reference_file "${local_reference_dir}/${reference_file}" )
  set( reference_diff_target "${target}_reference_diff" )
  set( local_reference_diff_target "${target}_bitid" )
  ectrans_test_fixture_name( target_fixture checksum_file ${target} )
  ectrans_local_reference_property_name( local_reference_property ${reference_name} )

  if( HAVE_REFERENCE_TESTS )
    ecbuild_add_test( TARGET ${reference_diff_target}
      COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/ectrans-diff.py"
      ARGS "${reference_file}" "${checksum_file}" -U 0 --max-output-lines 20 --overrideable-reference-root "${reference_dir}"
    )
    set_tests_properties( ${reference_diff_target} PROPERTIES
      DEPENDS ${target}
      FIXTURES_REQUIRED ${target_fixture}
      LABELS reference_test
      SKIP_RETURN_CODE 77 )
  endif()

  get_property( has_local_reference GLOBAL PROPERTY "${local_reference_property}_SOURCE_TARGET" SET )
  if( HAVE_BITID_TESTS AND has_local_reference )
    get_property( local_reference_source_target GLOBAL PROPERTY "${local_reference_property}_SOURCE_TARGET" )
    get_property( local_reference_source_fixture GLOBAL PROPERTY "${local_reference_property}_SOURCE_FIXTURE" )
    ecbuild_add_test( TARGET ${local_reference_diff_target}
      COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/ectrans-diff.py"
      ARGS "${local_reference_file}" "${checksum_file}" -U 0 --max-output-lines 20
    )
    set_tests_properties( ${local_reference_diff_target} PROPERTIES
      DEPENDS "${target};${local_reference_source_target}"
      FIXTURES_REQUIRED "${target_fixture};${local_reference_source_fixture}"
      LABELS "reference_test;local_reference"
      SKIP_RETURN_CODE 77 )
  endif()
endfunction()

# Apply common properties to an ectrans benchmark checksum producer.
# Use this on every benchmark test that may feed checksum comparison tests. It
# tags GPU tests and declares both the legacy aggregate checksum fixture and a
# per-test checksum fixture used by filtered reference-diff runs.
macro(ectrans_set_test_properties target)
  if( "${target}" MATCHES "gpu" )
    set_tests_properties(${target} PROPERTIES LABELS "gpu;Fortran")
  endif()
  ectrans_test_fixture_name( target_fixture checksum_file ${target} )
  set_tests_properties( ${target} PROPERTIES FIXTURES_SETUP "${target_fixture}" )
endmacro()
