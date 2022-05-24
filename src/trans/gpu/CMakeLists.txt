# (C) Copyright 2020- ECMWF.
# (C) Copyright 2022- NVIDIA.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                algor/*.F90
                                internal/*.F90
                                external/*.F90
                                internal/sharedmem.c
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
    ecbuild_add_library( 
      TARGET           trans_gpu_static_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_src}
      TYPE             STATIC
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
    ecbuild_add_library( 
      TARGET           trans_gpu_static_CA_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_src}
      TYPE             STATIC
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
    ectrans_target_fortran_module_directory(
      TARGET            trans_gpu_static_${prec}
      MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_gpu_static_${prec}
      INSTALL_DIRECTORY module/trans_gpu_static_${prec}
    )
    ectrans_target_fortran_module_directory(
      TARGET            trans_gpu_static_CA_${prec}
      MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_gpu_static_CA_${prec}
      INSTALL_DIRECTORY module/trans_gpu_static_CA_${prec}
    )

    if( HAVE_FFTW )
      target_link_libraries( trans_gpu_${prec} PRIVATE ${FFTW_LIBRARIES} )
      target_link_libraries( trans_gpu_static_${prec} PRIVATE ${FFTW_LIBRARIES} )
      target_link_libraries( trans_gpu_static_CA_${prec} PRIVATE ${FFTW_LIBRARIES} )
      target_include_directories( trans_gpu_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
      target_include_directories( trans_gpu_static_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
      target_include_directories( trans_gpu_static_CA_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
      target_compile_definitions( trans_gpu_${prec} PRIVATE WITH_FFTW )
    endif()

    if( HAVE_OMP )
      target_link_libraries( trans_gpu_${prec} PRIVATE OpenMP::OpenMP_Fortran )
      target_link_libraries( trans_gpu_static_${prec} PRIVATE OpenMP::OpenMP_Fortran )
      target_link_libraries( trans_gpu_static_CA_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    target_link_libraries( trans_gpu_${prec} PRIVATE OpenACC::OpenACC_Fortran )
    set_property( TARGET trans_gpu_${prec} PROPERTY CUDA_ARCHITECTURES 70 )
    target_compile_options( trans_gpu_${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-g -acc -Minfo=acc -gpu=cc70,lineinfo,deepcopy,fastmath,nordc> )
    target_compile_options( trans_gpu_static_${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-g -acc -Minfo=acc -gpu=cc70,lineinfo,deepcopy,fastmath,nordc> )
    target_compile_options( trans_gpu_static_CA_${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-g -acc -Minfo=acc -gpu=cc70,lineinfo,deepcopy,fastmath,nordc> )

    target_link_libraries( trans_gpu_static_${prec} PRIVATE OpenACC::OpenACC_Fortran )
    set_property( TARGET trans_gpu_static_${prec} PROPERTY CUDA_ARCHITECTURES 70 )

    target_link_libraries( trans_gpu_static_CA_${prec} PRIVATE OpenACC::OpenACC_Fortran )
    set_property( TARGET trans_gpu_static_CA_${prec} PROPERTY CUDA_ARCHITECTURES 70 )

    if( prec STREQUAL sp )
        target_compile_definitions( trans_gpu_${prec} PUBLIC TRANS_SINGLE PARKINDTRANS_SINGLE )
        target_compile_definitions( trans_gpu_static_${prec} PUBLIC TRANS_SINGLE PARKINDTRANS_SINGLE )
	target_compile_definitions( trans_gpu_static_CA_${prec} PUBLIC TRANS_SINGLE PARKINDTRANS_SINGLE )
    endif()
    target_compile_definitions( trans_gpu_static_CA_${prec} PUBLIC USE_CUDA_AWARE_MPI_FT )

  endif()
endforeach()

ecbuild_find_package( NAME NvidiaCutlass REQUIRED)

## precision-independent GPU library with CUDA kernels
ecbuild_add_library( TARGET gpu
  TYPE STATIC
  SOURCES
            algor/external/fourier/fft_wrapper.cu
            algor/external/gemm/gemm_wrapper.cu
            PRIVATE_INCLUDES  ${MPI_C_INCLUDE_PATH}
            PRIVATE_LIBS nvidia::cutlass::cutlass
)

## CUDA architecture 
set_property( TARGET gpu PROPERTY CUDA_ARCHITECTURES 70 ) 


## Install trans interface

file( GLOB trans_interface interface/* )
install(
  FILES        ${trans_interface}
  DESTINATION  include/ectrans
)
