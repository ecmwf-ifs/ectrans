# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Flags enabling OpenMP offload, applied only to the GPU targets through
# ectrans_target_omp_offload. OpenMP_Fortran_FLAGS is shared by every target linking
# OpenMP::OpenMP_Fortran, including the CPU libraries and other projects of a bundle (fiat),
# so it should keep host-only flags. FindOpenMP passes those flags to the compile step only
# for most compilers, so offload flags placed there leave non-ectrans executables unlinkable
# (nvfortran: undefined reference to `__acc_compiled').

# Default the flags per compiler unless set by the user. Cray needs none, as the
# craype-accel module enables offload.
if( GPU_OFFLOAD STREQUAL "OMP" AND NOT DEFINED ECTRANS_OMP_OFFLOAD_FLAGS )
  set( _ectrans_omp_offload_flags "" )
  if( CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
    set( _ectrans_omp_offload_flags "-mp=gpu" )
  elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang" )
    if( NOT CMAKE_HIP_ARCHITECTURES )
      ecbuild_critical( "OpenMP offload with ${CMAKE_Fortran_COMPILER_ID} needs CMAKE_HIP_ARCHITECTURES "
                        "or ECTRANS_OMP_OFFLOAD_FLAGS" )
    endif()
    # -fopenmp is repeated so the link, which gets only these flags, pulls in libomptarget.
    set( _ectrans_omp_offload_flags "-fopenmp" )
    foreach( _arch ${CMAKE_HIP_ARCHITECTURES} )
      list( APPEND _ectrans_omp_offload_flags "--offload-arch=${_arch}" )
    endforeach()
    list( JOIN _ectrans_omp_offload_flags " " _ectrans_omp_offload_flags )
  elseif( NOT CMAKE_Fortran_COMPILER_ID STREQUAL "Cray" )
    ecbuild_critical( "No default OpenMP offload flags for ${CMAKE_Fortran_COMPILER_ID}; "
                      "set ECTRANS_OMP_OFFLOAD_FLAGS" )
  endif()
  set( ECTRANS_OMP_OFFLOAD_FLAGS "${_ectrans_omp_offload_flags}" CACHE STRING
       "Flags enabling OpenMP offload for the ectrans GPU targets" )
  ecbuild_info( "ECTRANS_OMP_OFFLOAD_FLAGS: ${ECTRANS_OMP_OFFLOAD_FLAGS}" )
endif()

# Add ECTRANS_OMP_OFFLOAD_FLAGS to the Fortran compile of a GPU target, and to its link and that
# of everything linking it, as the final link must also enable offload.
function( ectrans_target_omp_offload target )
  if( GPU_OFFLOAD STREQUAL "OMP" AND ECTRANS_OMP_OFFLOAD_FLAGS )
    target_compile_options( ${target} PRIVATE
          $<$<COMPILE_LANGUAGE:Fortran>:SHELL:${ECTRANS_OMP_OFFLOAD_FLAGS}> )
    target_link_options( ${target} PUBLIC
          $<$<LINK_LANGUAGE:Fortran>:SHELL:${ECTRANS_OMP_OFFLOAD_FLAGS}>
          $<$<LINK_LANG_AND_ID:C,NVHPC>:SHELL:${ECTRANS_OMP_OFFLOAD_FLAGS}>
          $<$<LINK_LANG_AND_ID:CXX,NVHPC>:SHELL:${ECTRANS_OMP_OFFLOAD_FLAGS}> )
  endif()
endfunction()
