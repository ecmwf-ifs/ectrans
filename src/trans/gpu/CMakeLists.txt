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
  internal/tpm_stats.F90
  internal/tpm_hicfft.F90
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
  ecbuild_info("warn: IN_PLACE_FFT not defined for hipFFT")
elseif( HAVE_CUDA )
  set( GPU_RUNTIME "CUDA" )
  set( ECTRANS_GPU_HIP_LIBRARIES CUDA::cufft CUDA::cublas nvhpcwrapnvtx CUDA::cudart )
  list( APPEND trans_gpu_common_src
    algor/hicblas_gemm.cuda.cu
    algor/hicfft.cuda.cu
  )
  ecbuild_info("warn: IN_PLACE_FFT defined for cuFFT")
else()
  ecbuild_info("warn: HIP and CUDA not found")
endif()


set( GPU_LIBRARY_TYPE SHARED )
if( HAVE_GPU_STATIC )
    set( GPU_LIBRARY_TYPE STATIC )
endif()

ecbuild_add_library(
  TARGET               ectrans_gpu_common
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
                       $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                       $<${HAVE_CUTLASS}:nvidia::cutlass::cutlass>
  PRIVATE_DEFINITIONS  ${GPU_RUNTIME}GPU ${GPU_OFFLOAD}GPU
                       $<${HAVE_CUTLASS}:USE_CUTLASS>
                       $<${HAVE_CUTLASS_3XTF32}:USE_CUTLASS_3XTF32>
                       $<${HAVE_GPU_GRAPHS_GEMM}:USE_GRAPHS_GEMM>
                       $<${HAVE_GPU_GRAPHS_FFT}:USE_GRAPHS_FFT>
)

# The ecTrans libraries may be linked against an executable with the -cuda flag
# In that case we must link with -cuda here to ensure compatibility
target_link_options(ectrans_gpu_common PUBLIC
  $<$<LINK_LANG_AND_ID:Fortran,NVHPC>:-cuda>
)

ectrans_target_fortran_module_directory(
  TARGET            ectrans_gpu_common
  MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/ectrans
  INSTALL_DIRECTORY module/ectrans
)



function(generate_backend_sources)
  set (options)
  set (oneValueArgs BACKEND DESTINATION OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})
  file(MAKE_DIRECTORY ${destination}/algor)
  file(MAKE_DIRECTORY ${destination}/internal)
  file(MAKE_DIRECTORY ${destination}/external)

  ecbuild_list_add_pattern( LIST files
    GLOB
          internal/*.F90
          external/*.F90
    QUIET
  )
  list( APPEND files
          algor/seefmm_mix.F90
  )
  ecbuild_list_exclude_pattern( LIST files REGEX
          parkind_ectrans.F90
          tpm_stats.F90
          tpm_hicfft.F90
  )

  set(outfiles)
  foreach(file_i ${files})
    get_filename_component(outfile_name    ${file_i} NAME)
    get_filename_component(outfile_name_we ${file_i} NAME_WE)
    get_filename_component(outfile_ext     ${file_i} EXT)
    get_filename_component(outfile_dir     ${file_i} DIRECTORY)
    set(outfile "${destination}/${file_i}")
    ecbuild_debug("Generate ${outfile}")
    generate_file(BACKEND ${backend} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile})
    list(APPEND outfiles ${outfile})
  endforeach(file_i)
  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)


set( BUILD_INTERFACE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include/ectrans )

foreach( prec dp sp )
  if( HAVE_${prec} )

        set(GENERATED_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_gpu_${prec})

        generate_backend_includes(BACKEND gpu_${prec} TARGET ectrans_gpu_${prec}_includes DESTINATION ${BUILD_INTERFACE_INCLUDE_DIR} INCLUDE_DIRECTORY ${PROJECT_SOURCE_DIR}/src/trans/include )
        generate_backend_sources( BACKEND gpu_${prec} OUTPUT ectrans_gpu_${prec}_src  DESTINATION ${GENERATED_SOURCE_DIR})

        #if( NOT ${CMAKE_BUILD_TYPE_CAPS} STREQUAL DEBUG )
        set_source_files_properties( ${GENERATED_SOURCE_DIR}/internal/ftinv_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
        ecbuild_info("warn: special compile flags ftinv_mod.F90")
        set_source_files_properties( ${GENERATED_SOURCE_DIR}/internal/ftdir_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
        ecbuild_info("warn: special compile flags ftdir_mod.F90")
        #endif()


        ecbuild_add_library(
          TARGET               ectrans_gpu_${prec}
          TYPE                 ${GPU_LIBRARY_TYPE}
          SOURCES              ${ectrans_gpu_${prec}_src}
          LINKER_LANGUAGE      Fortran
          PUBLIC_INCLUDES      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
                               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                               $<INSTALL_INTERFACE:include/ectrans>
                               $<INSTALL_INTERFACE:include>
          PUBLIC_LIBS          ectrans_common ectrans_gpu_common ectrans_gpu_${prec}_includes
          PRIVATE_LIBS         ${ECTRANS_GPU_HIP_LIBRARIES}
                               $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                               $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                               $<${HAVE_MPI}:MPI::MPI_Fortran>
          PRIVATE_DEFINITIONS  ${GPU_RUNTIME}GPU ${GPU_OFFLOAD}GPU
                               $<${HAVE_CUTLASS}:USE_CUTLASS>
                               $<${HAVE_CUTLASS_3XTF32}:USE_CUTLASS_3XTF32>
                               $<${HAVE_GPU_GRAPHS_GEMM}:USE_GRAPHS_GEMM>
                               $<${HAVE_GPU_GRAPHS_FFT}:USE_GRAPHS_FFT>
                               $<${HAVE_GPU_AWARE_MPI}:USE_GPU_AWARE_MPI>
                               ECTRANS_HAVE_MPI=${ectrans_HAVE_MPI}
        )

        # The ecTrans libraries may be linked against an executable with the -cuda flag
        # In that case we must link with -cuda here to ensure compatibility
        target_link_options(ectrans_gpu_${prec} PUBLIC
          $<$<LINK_LANG_AND_ID:Fortran,NVHPC>:-cuda>
        )

        ectrans_target_fortran_module_directory(
          TARGET            ectrans_gpu_${prec}
          MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans
          INSTALL_DIRECTORY module/ectrans
        )

        if( prec STREQUAL sp )
          target_compile_definitions( ectrans_gpu_${prec} PRIVATE TRANS_SINGLE PARKINDTRANS_SINGLE )
        endif()

        # cuFFT can do in-place FFT, hipFFT cannot
        if( HAVE_CUDA )
          target_compile_definitions( ectrans_gpu_${prec} PRIVATE IN_PLACE_FFT )
        endif()

        if( HAVE_OMP AND CMAKE_Fortran_COMPILER_ID MATCHES Cray )
          # Propagate flags as link options for downstream targets. Only required for Cray
          target_link_options( ectrans_gpu_${prec} INTERFACE
                $<$<LINK_LANGUAGE:Fortran>:SHELL:${OpenMP_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:C,Cray>:SHELL:${OpenMP_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:CXX,Cray>:SHELL:${OpenMP_Fortran_FLAGS}> )
        endif()

        if( HAVE_ACC AND CMAKE_Fortran_COMPILER_ID MATCHES NVHPC )
          # Propagate flags as link options for downstream targets. Only required for NVHPC
          target_link_options( ectrans_gpu_${prec} INTERFACE
                $<$<LINK_LANGUAGE:Fortran>:SHELL:${OpenACC_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:C,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:CXX,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS}> )
        endif()

        # This interface library is for backward compatibility, and provides the older includes
        ecbuild_add_library( TARGET trans_gpu_${prec} TYPE INTERFACE )
        target_include_directories( trans_gpu_${prec} INTERFACE $<BUILD_INTERFACE:${BUILD_INTERFACE_INCLUDE_DIR}/trans_gpu_${prec}> )
        target_include_directories( trans_gpu_${prec} INTERFACE $<INSTALL_INTERFACE:include/ectrans/trans_gpu_${prec}> )
        target_link_libraries( trans_gpu_${prec} INTERFACE fiat ectrans_gpu_${prec} parkind_${prec})


        # ## Install trans_gpu_${prec} interface
        # file( GLOB trans_interface ${PROJECT_SOURCE_DIR}/src/trans/include/ectrans/* )
        # install(
        #   FILES        ${trans_interface}
        #   DESTINATION  include/ectrans/trans_gpu_${prec}
        # )
  endif()
endforeach()
