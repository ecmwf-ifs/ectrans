# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

function( ectrans_default_reference_dir out_var )
  set( _ectrans_reference_os "${CMAKE_SYSTEM_NAME}" )
  if( _ectrans_reference_os STREQUAL "Darwin" )
    set( _ectrans_reference_os "macos" )
  endif()

  set( _ectrans_reference_arch "${CMAKE_SYSTEM_PROCESSOR}" )
  if( NOT _ectrans_reference_arch )
    set( _ectrans_reference_arch "${CMAKE_HOST_SYSTEM_PROCESSOR}" )
  endif()

  set( _ectrans_reference_compiler "${CMAKE_Fortran_COMPILER_ID}" )
  if( NOT _ectrans_reference_compiler )
    set( _ectrans_reference_compiler "${ECBUILD_Fortran_COMPILER_ID}" )
  endif()
  if( NOT _ectrans_reference_compiler )
    set( _ectrans_reference_compiler "unknown" )
  endif()

  set( _ectrans_reference_build_type "${CMAKE_BUILD_TYPE}" )
  if( NOT _ectrans_reference_build_type )
    set( _ectrans_reference_build_type "default" )
  endif()

  set( _ectrans_reference_key
    "${_ectrans_reference_os}-${_ectrans_reference_arch}-${_ectrans_reference_compiler}-${_ectrans_reference_build_type}" )
  string( TOLOWER "${_ectrans_reference_key}" _ectrans_reference_key )
  string( REGEX REPLACE "[^a-z0-9]+" "-" _ectrans_reference_key "${_ectrans_reference_key}" )
  string( REGEX REPLACE "^-+|-+$" "" _ectrans_reference_key "${_ectrans_reference_key}" )

  set( ${out_var}
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/references/${_ectrans_reference_key}"
    PARENT_SCOPE )
endfunction()