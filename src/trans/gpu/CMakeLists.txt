# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                algor/*
                                internal/*
                                external/*
                                sharedmem/*
                          QUIET
                        )

ecbuild_list_exclude_pattern( LIST trans_src REGEX dilatation_mod.F90 )
if( NOT HAVE_FFTW )
  ecbuild_list_exclude_pattern( LIST trans_src REGEX tpm_fftw.F90 )
endif()


set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )
set( FFTW_dp FFTW::fftw3  )
set( FFTW_sp FFTW::fftw3f )

set( IFS_ACC_Fortran_LIBRARIES OpenACC::OpenACC_Fortran PARENT_SCOPE)

foreach( prec sp dp )
  if( HAVE_${prec} )

    ecbuild_add_library( 
      TARGET           trans_gpu_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_src}
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ectrans>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/algor/interface>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PRIVATE_INCLUDES  ${MPI_Fortran_INCLUDE_PATH}
      PUBLIC_LIBS      parkind_${prec}
                       fiat
      PRIVATE_LIBS     ${LAPACK_LIBRARIES}
    )
    ectrans_target_fortran_module_directory(
      TARGET            trans_gpu_${prec}
      MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_gpu_${prec}
      INSTALL_DIRECTORY module/trans_gpu_${prec}
    )
    if( HAVE_FFTW )
      target_link_libraries( trans_gpu_${prec} PRIVATE ${FFTW_LIBRARIES} )
      target_include_directories( trans_gpu_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
      target_compile_definitions( trans_gpu_${prec} PRIVATE WITH_FFTW )
    endif()
    if( HAVE_OMP )
      target_link_libraries( trans_gpu_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()
    target_link_libraries( trans_gpu_${prec} PRIVATE OpenACC::OpenACC_Fortran )
    set_property( TARGET trans_gpu_${prec} PROPERTY CUDA_ARCHITECTURES 70 )
    target_compile_options( trans_gpu_${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-g -acc -gpu=cc70,lineinfo,deepcopy,fastmath,nordc> )
  endif()
endforeach()

## precision-independent GPU library with CUDA kernels
ecbuild_add_library( TARGET gpu
  TYPE STATIC
  SOURCES
            algor/module/cublasSgemmBatched.cu
            algor/module/cublasDgemmBatched.cu
##            algor/module/cublasSTCgemmBatched.cu
            algor/module/IPC_Alltoall.cu
            algor/external/fourier/destroy_plan_fftc.cu
            algor/external/fourier/create_plan_fftc.cu
            algor/external/fourier/storage_fftc.cu
            algor/external/fourier/execute_plan_fftc.cu
            PRIVATE_INCLUDES  ${MPI_C_INCLUDE_PATH}
)

## CUDA architecture 
set_property( TARGET gpu PROPERTY CUDA_ARCHITECTURES 70 ) 


## Install trans interface

file( GLOB trans_interface interface/* )
install(
  FILES        ${trans_interface}
  DESTINATION  include/ectrans
)
