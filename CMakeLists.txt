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
	DEFAULT OFF
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
		    DEFAULT OFF
                    REQUIRED_PACKAGES "MKL QUIET" )

if( NOT HAVE_MKL )
    option( FFTW_ENABLE_MKL OFF )
endif()

if( ENABLE_GPU OR ECTRANS_ENABLE_GPU ) # must be explicitly turned ON as FEATURE GPU is OFF by default (see below)
        #include(CheckLanguage)
        #check_language(HIP)
        #if(CMAKE_HIP_COMPILER)
        #      enable_language(HIP)
        #else()
        #               ecbuild_critical("No HIP compiler found")
        #endif()
        #if(NOT DEFINED HIP_PATH)
        #if(NOT DEFINED ENV{HIP_PATH})
        #set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
        #else()
        #set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
        #endif()
        #endif()
        #list(APPEND CMAKE_PREFIX_PATH ${HIP_PATH} /opt/rocm)
        #set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
        set(CMAKE_MODULE_PATH $ENV{HIP_ROOT}/cmake ${CMAKE_MODULE_PATH})
        ecbuild_warn_var(CMAKE_MODULE_PATH)
        find_package(HIP)
        #if( HIP_FOUND )
        find_package(hipfft)
        find_package(hipblas)
        find_package(rocblas)
        find_package(rocfft)
        #endif()
  find_package(OpenACC COMPONENTS Fortran)
  find_package(MPI COMPONENTS Fortran CXX)

  if( NOT CMAKE_HIP_COMPILER)
    ecbuild_info("No HIP compiler found")
  endif()
  if(FALSE)
  if( NOT TARGET HIP::hipblas )
    ecbuild_info("No target HIP::hipblas")
  endif()
  if( NOT TARGET HIP::hipfft )
    ecbuild_info("No target HIP::hipfft")
  endif()
  endif()
  if( NOT TARGET OpenACC::OpenACC_Fortran )
    ecbuild_info("No target OpenACC:OpenACC_Fortran")
  endif()
  if( NOT TARGET MPI::MPI_Fortran )
    ecbuild_info("No target MPI::MPI_Fortran")
  endif()
  if( NOT TARGET MPI::MPI_CXX )
    ecbuild_info("No target MPI::MPI_CXX")
  endif()
endif()

ecbuild_add_option( FEATURE FFTW
                    DEFAULT ON
                    DESCRIPTION "Support for fftw"
                    REQUIRED_PACKAGES "FFTW COMPONENTS double ${single}" )

ecbuild_add_option( FEATURE TRANSI
                    DEFAULT ON
                    DESCRIPTION "Compile TransI C-interface to trans"
                    CONDITION HAVE_DOUBLE_PRECISION )

ecbuild_add_option( FEATURE GPU
                    DEFAULT OFF
                    DESCRIPTION "Compile GPU version of ectrans (Requires OpenACC or sufficient OpenMP offloading support and MPI and CUDA. For now only NVHPC supported.)"
                    CONDITION
                    HIP_FOUND
                         AND HAVE_MPI AND TARGET MPI::MPI_CXX AND TARGET MPI::MPI_Fortran
                         #     AND TARGET HIP::hipblas AND TARGET HIP::hipfft
                   )

ecbuild_add_option( FEATURE ACCGPU
                    DEFAULT ON
                    DESCRIPTION "Support for using GPUs with OpenACC"
                    CONDITION HAVE_GPU
                         AND TARGET OpenACC::OpenACC_Fortran)

ecbuild_add_option( FEATURE OMPGPU
                    DEFAULT ON
                    DESCRIPTION "Support for using GPUs with OpenMP offloading"
                    CONDITION HAVE_GPU
                         AND TARGET OpenMP::OpenMP_Fortran
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran" )

ecbuild_add_option( FEATURE GPU_AWARE_MPI
                    DEFAULT ON
                    CONDITION HAVE_GPU
                    DESCRIPTION "Enable CUDA-aware MPI")

ecbuild_add_option( FEATURE GPU_REDUCED_MEMORY
                    DEFAULT OFF
                    CONDITION HAVE_GPU
                    DESCRIPTION "Use (slower) code that reduces memory requirements on GPU" )

set( HAVE_acc ${HAVE_ACCGPU} )
set( HAVE_omp ${HAVE_OMPGPU} )

# Following is NVHPC compiler specific and should really be coming from external input
if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  if( HAVE_ACCGPU )
    set( acc_flags -acc -gpu=lineinfo,deepcopy,fastmath,nordc )
    set( acc_link_flags -acc )
    # Pass cmake command-line option "-DCMAKE_Fortran_FLAGS=-Minfo=acc" for diagnostics info
  endif()
  if( HAVE_OMPGPU )
    set( omp_flags -mp=gpu -gpu=lineinfo,fastmath -Minfo=mp )
    set( omp_link_flags -mp=gpu )
  endif()
endif()

if( HAVE_GPU )
    #enable_language( HIP )
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 70 80)
    endif()
endif()

#ectrans_find_lapack()

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
