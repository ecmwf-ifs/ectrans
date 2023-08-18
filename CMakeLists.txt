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

project( ectrans LANGUAGES C Fortran CXX )
include( ectrans_macros )

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
                    DEFAULT OFF
                    REQUIRED_PACKAGES "MKL QUIET" )

if( NOT HAVE_MKL )
    option( FFTW_ENABLE_MKL OFF )
endif()


ecbuild_add_option( FEATURE FFTW
                    DEFAULT ON
                    DESCRIPTION "Support for fftw"
                    REQUIRED_PACKAGES "FFTW COMPONENTS double ${single}" )

ecbuild_add_option( FEATURE CPU
                    DEFAULT ON
                    DESCRIPTION "Compile CPU version of ectrans"
                   )

ecbuild_add_option( FEATURE TRANSI
                    DEFAULT ON
                    DESCRIPTION "Compile TransI C-interface to trans"
                    CONDITION HAVE_DOUBLE_PRECISION AND HAVE_CPU )

set( HIP_REQUESTED OFF )
if( ECTRANS_ENABLE_GPU OR (NOT DEFINED ECTRANS_ENABLE_GPU AND ENABLE_GPU))
  set( HIP_REQUESTED ON )
endif()

if( HIP_REQUESTED )
  ectrans_find_hip() # sets "HAVE_HIP"
endif()           

ecbuild_add_option( FEATURE GPU
                    DEFAULT OFF
                    DESCRIPTION "Compile GPU version of ectrans (Requires OpenACC or sufficient OpenMP offloading support and MPI)"
                    CONDITION HAVE_HIP )

ecbuild_add_option( FEATURE ACC
                    DEFAULT ON
                    DESCRIPTION "Support for using GPUs with OpenACC"
                    REQUIRED_PACKAGES "OpenACC COMPONENTS Fortran"
                    CONDITION HAVE_GPU)
if( HAVE_GPU )
  if( HAVE_ACC )
    set( GPU_OFFLOAD "ACC" )
  elseif( HAVE_OMP )
    set( GPU_OFFLOAD "OMP" )
  else()
    ecbuild_error("Could not enable GPU as OMP or ACC were not enabled")
  endif()
endif()

ecbuild_add_option( FEATURE GPU_AWARE_MPI
                    DEFAULT OFF
                    CONDITION HAVE_GPU AND HAVE_MPI
                    DESCRIPTION "Enable CUDA-aware MPI" )

ecbuild_add_option( FEATURE GPU_REDUCED_MEMORY
                    DEFAULT OFF
                    CONDITION HAVE_GPU
                    DESCRIPTION "Use (slower) code that reduces memory requirements on GPU" )

if( BUILD_SHARED_LIBS )
  set( GPU_STATIC_DEFAULT OFF )
else()
  set( GPU_STATIC_DEFAULT ON )
endif()
ecbuild_add_option( FEATURE GPU_STATIC
                    DEFAULT ${GPU_STATIC_DEFAULT}
                    DESCRIPTION "Compile GPU library as static library")

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
