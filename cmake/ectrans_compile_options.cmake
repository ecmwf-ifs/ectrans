# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Flag to tell compiler that Fortran side has no program
# Needed if linking a C executable against some Fortran objects with some compilers
# Not needed for most
set( NO_FORTRAN_MAIN_FLAG "" )

if( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  ecbuild_add_fortran_flags("-qextname -qnobindcextname")
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  # gfortran 10 has become stricter with argument matching
  if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
    ecbuild_add_fortran_flags("-fallow-argument-mismatch")
  endif()
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  ecbuild_add_fortran_flags("-Mlarge_arrays")

  # should really be part of configuration, or ecbuild default?
  ecbuild_add_fortran_flags("-traceback"      BUILD DEBUG )
  ecbuild_add_fortran_flags("-fast"           BUILD RELEASE )
  ecbuild_add_fortran_flags("-gopt -fast"     BUILD RELWITHDEBINFO )

  set( NO_FORTRAN_MAIN_FLAG "-Mnomain")
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  ecbuild_add_fortran_flags("-hnomessage=878")  # A module named ... has already been directly or indirectly use associated into this scope
  ecbuild_add_fortran_flags("-hnomessage=867")  # Module ... has no public objects declared in the module, therefore nothing can be use associated from the module.
  ecbuild_add_fortran_flags("-M7256")           # An OpenMP parallel construct in a target region is limited to a single thread.
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM" )
  ecbuild_add_fortran_flags("-march=core-avx2 -no-fma" BUILD BIT)
  ecbuild_add_fortran_flags("-fp-model precise -fp-speculation=safe")
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  ecbuild_add_fortran_flags("-march=core-avx2 -no-fma" BUILD BIT)
  ecbuild_add_fortran_flags("-fast-transcendentals -fp-model precise -fp-speculation=safe")
  set( NO_FORTRAN_MAIN_FLAG "-nofor-main")
endif()

if( NOT DEFINED ECTRANS_HAVE_CONTIGUOUS_ISSUE )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel"  )
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS_EQUAL 19)
      set( ECTRANS_HAVE_CONTIGUOUS_ISSUE True )
    endif()
  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU"  )
    # GCC versions 9.2, 11.2, 12.2, 13.3, 14.2 are all known to have an issue with `contiguous`
    # Logic below is defensive and assumes future versions of gcc are likely to also have the issue
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 9 )
      set( ECTRANS_HAVE_CONTIGUOUS_ISSUE True )
    endif()
  endif()
endif()

macro( ectrans_add_compile_options )
  set( options NOFAIL )
  set( single_value_args FLAGS)
  set( multi_value_args SOURCES )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )
  if(_PAR_UNPARSED_ARGUMENTS)
    ecbuild_critical("Unknown keywords given to ectrans_add_compile_flags(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
  endif()
  if(NOT _PAR_SOURCES)
    ecbuild_critical("SOURCES keyword missing to ectrans_add_compile_flags()")
  endif()
  if(NOT _PAR_FLAGS)
    ecbuild_critical("FLAGS keyword missing to ectrans_add_compile_flags()")
  endif()
  foreach( _file ${_PAR_SOURCES} )
    ecbuild_warn("Adding custom compile flags for file ${_file} : [${_PAR_FLAGS}]")
    if( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_file} AND NOT _PAR_NOFAIL)
        ecbuild_error("${_file} does not exist")
    endif()
    set_source_files_properties( ${_file} PROPERTIES COMPILE_FLAGS "${_PAR_FLAGS}" )
  endforeach()
endmacro()

