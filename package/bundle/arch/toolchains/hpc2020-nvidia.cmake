####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI ON )
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_C_COMPILER nvc)
set(CMAKE_CXX_COMPILER nvc++)
set(CMAKE_Fortran_COMPILER nvfortran)

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-mp -mp=bind,allcores,numa" )
set( OpenMP_CXX_FLAGS           "-mp -mp=bind,allcores,numa" )
if(ENABLE_OMP_OFFLOAD)
  set( OpenMP_Fortran_FLAGS       "-mp -mp=gpu,bind,allcores,numa -gpu=cc80,lineinfo,fastmath,rdc" CACHE STRING "" FORCE)
else()
  set( OpenMP_Fortran_FLAGS       "-mp -mp=bind,allcores,numa" CACHE STRING "" FORCE)
endif()


####################################################################
# OpenAcc FLAGS
####################################################################

set( OpenACC_Fortran_FLAGS "-acc=gpu -gpu=lineinfo,fastmath,rdc" )

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "-fpic")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mframe")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mstack_arrays")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mrecursive")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kieee")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mdaz")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -I/usr/local/apps/hpcx-openmpi/2.14.0-cuda/NVIDIA/22.11/ec-hpcx-ompi/lib")

set( ECBUILD_Fortran_FLAGS_BIT "-O2 -gopt" )

set( ECBUILD_C_FLAGS "-O2 -gopt -traceback" )

set( ECBUILD_CXX_FLAGS "-O2 -gopt" )
