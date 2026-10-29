# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# Compiler
####################################################################

set(CMAKE_C_COMPILER nvc)
set(CMAKE_CXX_COMPILER nvc++)
set(CMAKE_Fortran_COMPILER nvfortran)

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-mp -mp=bind,allcores,numa" )
set( OpenMP_CXX_FLAGS           "-mp -mp=bind,allcores,numa" )
set( OpenMP_Fortran_FLAGS       "-mp -mp=bind,allcores,numa" )

# Added by ectrans to its GPU targets only, when built with OpenMP offload
set( ECTRANS_OMP_OFFLOAD_FLAGS  "-mp=gpu -gpu=cc80,lineinfo,fastmath,rdc" )

####################################################################
# OpenACC FLAGS
####################################################################

set( OpenACC_Fortran_FLAGS "-acc=gpu -gpu=cc80,lineinfo,fastmath,rdc" )
set( ECTRANS_ACC_UNIFIED_FLAGS "-gpu=mem:unified:nomanagedalloc" )
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

####################################################################
# Default features
####################################################################

set( ENABLE_GPU ON )
if( NOT DEFINED ENABLE_ACC )
  set( ENABLE_ACC ON )
endif()
set( ENABLE_GPU_STATIC ON )
set( ENABLE_FIELD_API ON )
