# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Preprocess module file containing version information
configure_file( internal/ectrans_version_mod.F90.in internal/ectrans_version_mod.F90 )

## Apply workarounds for some known compilers

if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")

  # Compile setup_trans with pinned memory to improve data movement performance.
  ectrans_add_compile_options(
      SOURCES external/setup_trans.F90
      #FLAGS   "-gpu=pinned,deepcopy,fastmath,nordc")
      FLAGS   "-gpu=pinned,fastmath")
  # TODO: check if it is sufficient to only set "-gpu=pinned" which appends rather than overwrites

endif()

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                sharedmem/*
                                algor/*
                                internal/*
                                external/*
                                ${CMAKE_CURRENT_BINARY_DIR}/internal/ectrans_version_mod.F90
                          QUIET
                        )

ecbuild_list_exclude_pattern( LIST trans_src            REGEX dilatation_mod.F90 )

#if( NOT ${CMAKE_BUILD_TYPE_CAPS} STREQUAL DEBUG )
  set_source_files_properties( internal/ftinv_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
  ecbuild_info("warn: special compile flags ftinv_mod.F90")
  set_source_files_properties( internal/ftdir_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
  ecbuild_info("warn: special compile flags ftdir_mod.F90")
#endif()

# Filter source list according to available GPU runtime
if( HAVE_HIP )
  set( GPU_RUNTIME "HIP" )
  ectrans_declare_hip_sources( SOURCES_GLOB
    sharedmem/*.hip.cpp
    algor/*.hip.cpp
    internal/*.hip.cpp
    external/*.hip.cpp
  )
  ecbuild_list_exclude_pattern( LIST trans_src  REGEX \.cu$ )
  ecbuild_list_exclude_pattern( LIST trans_src  REGEX cuda_device_mod.F90 )
elseif( HAVE_CUDA )
  set( GPU_RUNTIME "CUDA" )
  set( ECTRANS_GPU_HIP_LIBRARIES CUDA::cufft CUDA::cublas nvhpcwrapnvtx CUDA::cudart )
  ecbuild_list_exclude_pattern( LIST trans_src  REGEX \.hip\.cpp )
  ecbuild_list_exclude_pattern( LIST trans_src  REGEX hip_device_mod.F90 )
else()
  ecbuild_info("warn: HIP and CUDA not found")
endif()

foreach( prec dp sp )
  if( HAVE_${prec} )
        set( GPU_LIBRARY_TYPE SHARED )
        if( HAVE_GPU_STATIC )
            set( GPU_LIBRARY_TYPE STATIC )
        endif()

        ecbuild_add_library(
          TARGET               trans_gpu_${prec}
          TYPE                 ${GPU_LIBRARY_TYPE}
          SOURCES              ${trans_src}
          LINKER_LANGUAGE      Fortran
          PUBLIC_INCLUDES      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
                               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/algor/interface>
                               $<INSTALL_INTERFACE:include/ectrans/trans_gpu_${prec}>
                               $<INSTALL_INTERFACE:include/ectrans>
                               $<INSTALL_INTERFACE:include>
          PUBLIC_LIBS          parkind_${prec}
                               fiat
          PRIVATE_LIBS         ${ECTRANS_GPU_HIP_LIBRARIES}
                               ${LAPACK_LIBRARIES} # we still have symbols in some files
                               $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                               $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                               MPI::MPI_Fortran MPI::MPI_CXX
                               $<${HAVE_CUTLASS}:nvidia::cutlass::cutlass>
          PRIVATE_DEFINITIONS  ${GPU_OFFLOAD}GPU
                               ${GPU_RUNTIME}GPU
                               $<${HAVE_GPU_AWARE_MPI}:USE_GPU_AWARE_MPI>
                               $<${HAVE_CUTLASS}:USE_CUTLASS>
                               $<${HAVE_CUTLASS_3XTF32}:USE_CUTLASS_3XTF32>
                               $<${HAVE_GPU_GRAPHS_GEMM}:USE_GRAPHS_GEMM>
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
