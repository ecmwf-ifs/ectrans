# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )
find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

project( ectrans LANGUAGES C Fortran )
include( ectrans_macros )

ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

### Find (optional) dependencies 

ecbuild_find_package( NAME fiat REQUIRED )

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for MPI distributed memory parallelism"
                    CONDITION fiat_HAVE_MPI )

ecbuild_add_option( FEATURE OMP
                    DEFAULT ON
                    DESCRIPTION "Support for OpenMP shared memory parallelism"
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran" )

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

ecbuild_add_option( FEATURE FFTW
                    DEFAULT ON
                    DESCRIPTION "Support for fftw"
                    REQUIRED_PACKAGES "FFTW COMPONENTS double ${single}" )

ecbuild_add_option( FEATURE TRANSI
                    DEFAULT ON
                    DESCRIPTION "Compile TransI C-interface to trans"
                    CONDITION HAVE_DOUBLE_PRECISION )

# Following is NVHPC compiler specific and should really be coming from external input
if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  set( OpenACC_Fortran_FLAGS "-acc -gpu=lineinfo,deepcopy,fastmath,nordc" )
  # Pass cmake command-line option "-DCMAKE_Fortran_FLAGS=-Minfo=acc" for diagnostics info
endif()

if( ENABLE_GPU OR ECTRANS_ENABLE_GPU ) # must be explicitly turned ON as FEATURE GPU is OFF by default (see below)
  check_language(CUDA)
  find_package(CUDAToolkit)
  find_package(OpenACC COMPONENTS Fortran)
  find_package(MPI COMPONENTS Fortran C)

  if( NOT CMAKE_CUDA_COMPILER)
    ecbuild_info("No CUDA compiler found")
  endif()
  if( NOT TARGET CUDA::cublas )
    ecbuild_info("No target CUDA::cublas")
  endif()
  if( NOT TARGET CUDA::cufft )
    ecbuild_info("No target CUDA::cufft")
  endif()
  if( NOT TARGET OpenACC::OpenACC_Fortran )
    ecbuild_info("No target OpenACC:OpenACC_Fortran")
  endif()
  if( NOT TARGET MPI::MPI_Fortran )
    ecbuild_info("No target MPI::MPI_Fortran")
  endif()
  if( NOT TARGET MPI::MPI_C )
    ecbuild_info("No target MPI::MPI_C")
  endif()
endif()

ecbuild_add_option( FEATURE GPU
                    DEFAULT OFF
                    DESCRIPTION "Compile GPU version of ectrans (Requires OpenACC and MPI and CUDA. For now only NVHPC supported.)"
                    CONDITION CMAKE_CUDA_COMPILER
                         AND CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC"
                         AND TARGET OpenACC::OpenACC_Fortran
                         AND HAVE_MPI AND TARGET MPI::MPI_C AND TARGET MPI::MPI_Fortran
                         AND TARGET CUDA::cublas AND TARGET CUDA::cufft)

ecbuild_add_option( FEATURE GPU_AWARE_MPI
                    DEFAULT ON
                    CONDITION HAVE_GPU
                    DESCRIPTION "Enable CUDA-aware MPI")

ecbuild_add_option( FEATURE GPU_REDUCED_MEMORY
                    DEFAULT OFF
                    CONDITION HAVE_GPU
                    DESCRIPTION "Use (slower) code that reduces memory requirements on GPU" )

if( HAVE_GPU )
    enable_language( CUDA )
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 70 80)
    endif()
endif()

ectrans_find_lapack()

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

