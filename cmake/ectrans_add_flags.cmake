# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ectrans_add_flags flags )
  cmake_parse_arguments( _PAR "" "BUILD;LANG" "" ${ARGN} )
  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER )
  string(TOUPPER "${_PAR_BUILD}" _PAR_BUILD_UPPER)
  if( _PAR_BUILD )
    if (CMAKE_BUILD_TYPE_UPPER STREQUAL _PAR_BUILD_UPPER)
      ecbuild_info( "Adding ${flags} to ${_PAR_LANG} flags for build type ${CMAKE_BUILD_TYPE}" )
    endif()
    set( FLAGS_VAR_NAME "${PNAME}_${_PAR_LANG}_FLAGS_${_PAR_BUILD_UPPER}" )
  else()
    ecbuild_info( "Adding ${flags} to ${_PAR_LANG} flags for all build types" )
    set( FLAGS_VAR_NAME "${PNAME}_${_PAR_LANG}_FLAGS" )
  endif()
  if ( NOT ${FLAGS_VAR_NAME} )
    set( ${FLAGS_VAR_NAME} "${flags}" )
  else()
    set( ${FLAGS_VAR_NAME} "${${FLAGS_VAR_NAME}} ${flags}" )
  endif()
endmacro()

macro (ectrans_add_fortran_flags flags )
  ectrans_add_flags( ${flags} LANG Fortran ${ARGN} )
endmacro()

macro (ectrans_add_c_flags flags )
  ectrans_add_flags( ${flags} LANG C ${ARGN} )
endmacro()

macro (ectrans_add_cxx_flags flags )
  ectrans_add_flags( ${flags} LANG CXX ${ARGN} )
endmacro()