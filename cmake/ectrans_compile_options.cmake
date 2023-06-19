# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


if( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  ecbuild_add_fortran_flags("-qextname -qnobindcextname")
endif()

# gfortran 10 has become stricter with argument matching
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU"
    AND NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
  ecbuild_add_fortran_flags("-fallow-argument-mismatch")
endif()

if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  ecbuild_add_fortran_flags("-Mlarge_arrays")

  # should really be part of configuration, or ecbuild default?
  ecbuild_add_fortran_flags("-traceback"      BUILD DEBUG )
  ecbuild_add_fortran_flags("-fast"           BUILD RELEASE )
  ecbuild_add_fortran_flags("-gopt -fast"     BUILD RELWITHDEBINFO )
endif()

if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  ecbuild_add_fortran_flags("-hnomessage=878")  # A module named ... has already been directly or indirectly use associated into this scope
  ecbuild_add_fortran_flags("-hnomessage=867")  # Module ... has no public objects declared in the module, therefore nothing can be use associated from the module.
  ecbuild_add_fortran_flags("-M7256")           # An OpenMP parallel construct in a target region is limited to a single thread.
endif()


############ !!!!! These should not be set within the project !!!!!
# Following is NVHPC compiler specific and should really be coming from external input
# if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
#   if( HAVE_ACCGPU )
#     set( acc_flags -acc -gpu=lineinfo,deepcopy,fastmath,nordc )
#     set( acc_link_flags -acc )
#     # Pass cmake command-line option "-DCMAKE_Fortran_FLAGS=-Minfo=acc" for diagnostics info
#   endif()
#   if( HAVE_OMPGPU )
#     set( omp_flags -mp=gpu -gpu=lineinfo,fastmath -Minfo=mp )
#     set( omp_link_flags -mp=gpu )
#   endif()
# endif()
#
# Following is Cray compiler specific and should really be coming from external input
# if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
#   if( HAVE_ACCGPU )
#     set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -hacc_model=auto_async_kernel:no_fast_addr:deep_copy")
#     set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -hacc -hacc_model=auto_async_kernel:no_fast_addr:deep_copy")
#   endif()
#   if( HAVE_OMPGPU )
#     set( omp_flags -fopenmp )
#     set( omp_link_flags -fopenmp )
#   endif()
# endif()
############




macro( ectrans_add_compile_options )
  set( options  )
  set( single_value_args FLAGS )
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
    if( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_file} )
        ecbuild_error("${_file} does not exist")
    endif()
    set_source_files_properties( ${_file} PROPERTIES COMPILE_FLAGS "${_PAR_FLAGS}" )
  endforeach()
endmacro()

