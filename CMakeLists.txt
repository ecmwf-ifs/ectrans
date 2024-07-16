# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.18 FATAL_ERROR )
# CMake 3.17 adds INTERFACE link options which get propagated to the link stage,
#   even if the target is linked in privately
# CMake 3.18 allows "LINK_LANG_AND_ID" generator expression.

find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

project( ectrans LANGUAGES C CXX Fortran )
include( ectrans_macros )

set(CMAKE_CXX_STANDARD 17)

option(IFS_SHARED_LIBS "Build shared libraries (overrides BUILD_SHARED_LIBS)" ON)
set(BUILD_SHARED_LIBS ${IFS_SHARED_LIBS})

ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

### Find (optional) dependencies

ecbuild_find_package( NAME fiat REQUIRED )

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for MPI distributed memory parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS Fortran CXX"
                    CONDITION fiat_HAVE_MPI )

ecbuild_add_option( FEATURE OMP
                    DEFAULT ON
                    DESCRIPTION "Support for OpenMP shared memory parallelism"
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran" )

if( ${CMAKE_VERSION} VERSION_LESS "3.25" AND (NOT DEFINED ENABLE_ACC OR ENABLE_ACC ) )
  # See https://gitlab.kitware.com/cmake/cmake/-/issues/23691, fixed in CMake 3.25
  # (TL;DR: FindOpenACC sets OpenACC_<LANG>_FOUND correctly but does not set
  #  OpenACC_FOUND unless all three C, CXX, and Fortran have been found - even if
  #  only one language has been requested via COMPONENTS)
  find_package( OpenACC COMPONENTS Fortran )
  if( OpenACC_Fortran_FOUND )
    set( OpenACC_FOUND ON )
  endif()
endif()

ecbuild_add_option( FEATURE ACC
                    DEFAULT OFF
                    DESCRIPTION "Support for using GPUs with OpenACC"
                    REQUIRED_PACKAGES "OpenACC COMPONENTS Fortran" )

ecbuild_add_option( FEATURE DOUBLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Support for Double Precision" )

ecbuild_add_option( FEATURE SINGLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Support for Single Precision" )

if( HAVE_SINGLE_PRECISION )
  set( single "single" )
endif()
set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

ecbuild_add_option( FEATURE MKL
                    DESCRIPTION "Use MKL for BLAS and/or FFTW"
                    DEFAULT ON
                    REQUIRED_PACKAGES "MKL QUIET" )

if( NOT HAVE_MKL )
    option( FFTW_ENABLE_MKL OFF )
endif()

ecbuild_add_option( FEATURE CPU
                    DEFAULT ON
                    DESCRIPTION "Compile CPU version of ectrans"
                   )

if( HAVE_CPU )
  ecbuild_find_package( NAME FFTW REQUIRED COMPONENTS double ${single} )
endif()

ecbuild_add_option( FEATURE TRANSI
                    DEFAULT ON
                    DESCRIPTION "Compile TransI C-interface to trans"
                    CONDITION HAVE_DOUBLE_PRECISION AND HAVE_CPU )

# Search for available GPU runtimes, searching for CUDA first and, if not found,
# attempt to find HIP
if( ECTRANS_ENABLE_GPU OR (NOT DEFINED ECTRANS_ENABLE_GPU AND ENABLE_GPU))
  set(HAVE_CUDA 0)
  set(HAVE_HIP 0)
  ectrans_find_cuda() # sets "HAVE_CUDA"
  if( NOT HAVE_CUDA )
    if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
        ectrans_find_hip() # sets "HAVE_HIP"
    else()
        message("Cannot enable for HIP language, require CMake >= 3.24")
    endif()
  endif()
endif()

ecbuild_add_option( FEATURE GPU
                    DEFAULT OFF
                    DESCRIPTION "Compile GPU version of ectrans (Requires OpenACC or sufficient OpenMP offloading support and MPI)"
                    CONDITION (HAVE_HIP OR HAVE_CUDA) AND (HAVE_ACC OR HAVE_OMP) AND HAVE_MPI )

if( HAVE_GPU )
  if( HAVE_ACC )
    set( GPU_OFFLOAD "ACC" )
  elseif( HAVE_OMP )
    set( GPU_OFFLOAD "OMP" )
  else()
    ecbuild_error("Could not enable GPU as OMP or ACC were not enabled")
  endif()
endif()

ecbuild_add_option( FEATURE CUTLASS
                    DEFAULT OFF
                    CONDITION HAVE_GPU AND HAVE_CUDA AND CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC"
                    DESCRIPTION "Support for Cutlass BLAS operations"
                    REQUIRED_PACKAGES "NvidiaCutlass VERSION 2.11" )

# following also needs cuda arch sm80 to be effective
ecbuild_add_option( FEATURE CUTLASS_3XTF32
                    DEFAULT ON
                    CONDITION HAVE_SINGLE_PRECISION AND HAVE_CUTLASS
                    DESCRIPTION "Support for 3xTF32 with Cutlass (>= 2.8) and CUDA_ARCHITECTURES >= 80" )

ecbuild_add_option( FEATURE GPU_AWARE_MPI
                    DEFAULT ON
                    CONDITION HAVE_GPU
                    REQUIRED_PACKAGES "MPI COMPONENTS CXX Fortran"
                    DESCRIPTION "Enable CUDA-aware MPI" )

ecbuild_add_option( FEATURE GPU_GRAPHS_GEMM
                    DEFAULT ON
                    CONDITION HAVE_GPU
                    DESCRIPTION "Enable graph-based optimisation of Legendre transform GEMM kernel" )

if( BUILD_SHARED_LIBS )
  set( GPU_STATIC_DEFAULT OFF )
else()
  set( GPU_STATIC_DEFAULT ON )
endif()
ecbuild_add_option( FEATURE GPU_STATIC
                    DEFAULT ${GPU_STATIC_DEFAULT}
                    DESCRIPTION "Compile GPU library as static library")

ectrans_find_lapack()

ecbuild_add_option( FEATURE TESTS
                    DEFAULT ON
                    DESCRIPTION "Enable unit testing"
                    REQUIRED_PACKAGES "MPI COMPONENTS Fortran"
                    CONDITION HAVE_CPU )

### Add sources and tests
include( ectrans_compile_options )
add_subdirectory( src )
add_subdirectory( tests )

### Export
if( BUILD_SHARED_LIBS )
  set( PACKAGE_REQUIRES_PRIVATE_DEPENDENCIES 0 )
else()
  set( PACKAGE_REQUIRES_PRIVATE_DEPENDENCIES 1 )
endif()

ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()
