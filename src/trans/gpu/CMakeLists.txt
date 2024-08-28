# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Assemble sources

list( APPEND trans_gpu_common_src
  algor/ext_acc.F90
  algor/c_hipmemgetinfo.cpp
  algor/buffered_allocator_mod.F90
  algor/device_mod.F90
  algor/growing_allocator_mod.F90
  algor/hicblas_mod.F90
)
if( HAVE_HIP )
  set( GPU_RUNTIME "HIP" )
  ectrans_declare_hip_sources( SOURCES_GLOB
    algor/*.hip.cpp
  )
  list( APPEND trans_gpu_common_src
    algor/hicblas_gemm.hip.cpp
    algor/hicfft.hip.cpp
  )
elseif( HAVE_CUDA )
  set( GPU_RUNTIME "CUDA" )
  set( ECTRANS_GPU_HIP_LIBRARIES CUDA::cufft CUDA::cublas nvhpcwrapnvtx CUDA::cudart )
  list( APPEND trans_gpu_common_src
    algor/hicblas_gemm.cuda.cu
    algor/hicfft.cuda.cu
  )
else()
  ecbuild_info("warn: HIP and CUDA not found")
endif()


set( GPU_LIBRARY_TYPE SHARED )
if( HAVE_GPU_STATIC )
    set( GPU_LIBRARY_TYPE STATIC )
endif()

ecbuild_add_library(
  TARGET               trans_gpu_common
  TYPE                 ${GPU_LIBRARY_TYPE}
  SOURCES              ${trans_gpu_common_src}
  LINKER_LANGUAGE      Fortran
  PUBLIC_INCLUDES      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
                       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
  PUBLIC_LIBS          fiat ectrans_common
  PRIVATE_LIBS         ${ECTRANS_GPU_HIP_LIBRARIES}
                       $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                       $<${HAVE_CUTLASS}:nvidia::cutlass::cutlass>
  PRIVATE_DEFINITIONS  ${GPU_RUNTIME}GPU
                       $<${HAVE_CUTLASS}:USE_CUTLASS>
                       $<${HAVE_CUTLASS_3XTF32}:USE_CUTLASS_3XTF32>
                       $<${HAVE_GPU_GRAPHS_GEMM}:USE_GRAPHS_GEMM>
)

ectrans_target_fortran_module_directory(
  TARGET            trans_gpu_common
  MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/ectrans
  INSTALL_DIRECTORY module/ectrans
)



ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                internal/*
                                external/*
                          QUIET
                        )

#if( NOT ${CMAKE_BUILD_TYPE_CAPS} STREQUAL DEBUG )
  set_source_files_properties( internal/ftinv_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
  ecbuild_info("warn: special compile flags ftinv_mod.F90")
  set_source_files_properties( internal/ftdir_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
  ecbuild_info("warn: special compile flags ftdir_mod.F90")
#endif()

foreach( prec dp sp )
  if( HAVE_${prec} )
        ecbuild_add_library(
          TARGET               trans_gpu_${prec}
          TYPE                 ${GPU_LIBRARY_TYPE}
          SOURCES              ${trans_src}
                               algor/seefmm_mix.F90
                               algor/wts500_mod.F90
          LINKER_LANGUAGE      Fortran
          PUBLIC_INCLUDES      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
                               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                               $<INSTALL_INTERFACE:include/ectrans/trans_gpu_${prec}>
                               $<INSTALL_INTERFACE:include/ectrans>
                               $<INSTALL_INTERFACE:include>
          PUBLIC_LIBS          parkind_${prec}
                               fiat
                               ectrans_common
          PRIVATE_LIBS         trans_gpu_common
                               ${ECTRANS_GPU_HIP_LIBRARIES}
                               $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                               $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                               $<${HAVE_MPI}:MPI::MPI_Fortran>
          PRIVATE_DEFINITIONS  ${GPU_OFFLOAD}GPU
                               ECTRANS_HAVE_MPI=${ectrans_HAVE_MPI}
                               $<${HAVE_GPU_AWARE_MPI}:USE_GPU_AWARE_MPI>
        )

        ectrans_target_fortran_module_directory(
          TARGET            trans_gpu_${prec}
          MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_gpu_${prec}
          INSTALL_DIRECTORY module/trans_gpu_${prec}
        )

        if( prec STREQUAL sp )
          target_compile_definitions( trans_gpu_${prec} PRIVATE TRANS_SINGLE PARKINDTRANS_SINGLE )
        endif()

        if( HAVE_OMP AND CMAKE_Fortran_COMPILER_ID MATCHES Cray )
          # Propagate flags as link options for downstream targets. Only required for Cray
          target_link_options( trans_gpu_${prec} INTERFACE
                $<$<LINK_LANGUAGE:Fortran>:SHELL:${OpenMP_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:C,Cray>:SHELL:${OpenMP_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:CXX,Cray>:SHELL:${OpenMP_Fortran_FLAGS}> )
        endif()

        if( HAVE_ACC AND CMAKE_Fortran_COMPILER_ID MATCHES NVHPC )
          # Propagate flags as link options for downstream targets. Only required for NVHPC
          target_link_options( trans_gpu_${prec} INTERFACE
                $<$<LINK_LANGUAGE:Fortran>:SHELL:${OpenACC_Fortran_FLAGS} -gpu=pinned>
                $<$<LINK_LANG_AND_ID:C,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS} -gpu=pinned>
                $<$<LINK_LANG_AND_ID:CXX,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS} -gpu=pinned> )
        endif()

        ## Install trans_gpu_${prec} interface
        file( GLOB trans_interface ${PROJECT_SOURCE_DIR}/src/trans/include/ectrans/* )
        install(
          FILES        ${trans_interface}
          DESTINATION  include/ectrans/trans_gpu_${prec}
        )
  endif()
endforeach()
