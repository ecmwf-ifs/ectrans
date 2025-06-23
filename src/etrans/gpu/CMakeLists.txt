set( GPU_LIBRARY_TYPE SHARED )
if( HAVE_GPU_STATIC )
    set( GPU_LIBRARY_TYPE STATIC )
endif()

function(generate_backend_sources)
  set (options)
  set (oneValueArgs BACKEND DESTINATION OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})
  file(MAKE_DIRECTORY ${destination}/internal)
  file(MAKE_DIRECTORY ${destination}/external)
  file(MAKE_DIRECTORY ${destination}/biper/external)
  file(MAKE_DIRECTORY ${destination}/biper/internal)

  ecbuild_list_add_pattern( LIST files
    GLOB
          internal/*.F90
          external/*.F90
          biper/internal/*.F90
          biper/external/*.F90
    QUIET
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

        set(GENERATED_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_lam_gpu_${prec})

        generate_backend_includes(BACKEND gpu_${prec} TARGET ectrans_lam_gpu_${prec}_includes DESTINATION ${BUILD_INTERFACE_INCLUDE_DIR} INCLUDE_DIRECTORY ${PROJECT_SOURCE_DIR}/src/etrans/include )
        generate_backend_sources( BACKEND gpu_${prec} OUTPUT ectrans_lam_gpu_${prec}_src  DESTINATION ${GENERATED_SOURCE_DIR})


        # set custom compilation flags here: keeping as placeholder
        #if( NOT ${CMAKE_BUILD_TYPE_CAPS} STREQUAL DEBUG )
        #set_source_files_properties( ${GENERATED_SOURCE_DIR}/internal/ftinv_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
        #ecbuild_info("warn: special compile flags ftinv_mod.F90")
        #set_source_files_properties( ${GENERATED_SOURCE_DIR}/internal/ftdir_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
        #ecbuild_info("warn: special compile flags ftdir_mod.F90")
        #endif()

        ecbuild_add_library(
          TARGET               ectrans_lam_gpu_${prec}
          TYPE                 ${GPU_LIBRARY_TYPE}
          SOURCES              ${ectrans_lam_gpu_${prec}_src}
          LINKER_LANGUAGE      Fortran
          PUBLIC_INCLUDES      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/etrans/include>
                               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/etrans/include/etrans>
                               $<INSTALL_INTERFACE:include/ectrans>
                               $<INSTALL_INTERFACE:include>
          PUBLIC_LIBS          fiat ectrans_common ectrans_gpu_common ectrans_gpu_${prec} ectrans_gpu_${prec}_includes ectrans_lam_common ectrans_lam_gpu_${prec}_includes
          PRIVATE_LIBS         ${ECTRANS_GPU_HIP_LIBRARIES}
                               $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                               $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                               $<${HAVE_MPI}:MPI::MPI_Fortran>
          PRIVATE_DEFINITIONS  ${GPU_RUNTIME}GPU ${GPU_OFFLOAD}GPU
                               #$<${HAVE_CUTLASS}:USE_CUTLASS>   # not relevant for LAM
                               #$<${HAVE_CUTLASS_3XTF32}:USE_CUTLASS_3XTF32>   # not relevant for LAM
                               #$<${HAVE_GPU_GRAPHS_GEMM}:USE_GRAPHS_GEMM>   # not relevant for LAM
                               $<${HAVE_GPU_GRAPHS_FFT}:USE_GRAPHS_FFT>
                               #$<${HAVE_GPU_AWARE_MPI}:USE_GPU_AWARE_MPI>   # not relevant for LAM
                               #ECTRANS_HAVE_MPI=${ectrans_HAVE_MPI}   # not relevant for LAM
        )

        ectrans_target_fortran_module_directory(
          TARGET            ectrans_lam_gpu_${prec}
          MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans
          INSTALL_DIRECTORY module/ectrans
        )

        if( prec STREQUAL sp )
          target_compile_definitions( ectrans_lam_gpu_${prec} PRIVATE TRANS_SINGLE PARKINDTRANS_SINGLE )
        endif()

        # cuFFT can do in-place FFT, hipFFT cannot
        if( HAVE_CUDA )
          target_compile_definitions( ectrans_lam_gpu_${prec} PRIVATE IN_PLACE_FFT )
        endif()

        if( HAVE_OMP AND CMAKE_Fortran_COMPILER_ID MATCHES Cray )
          # Propagate flags as link options for downstream targets. Only required for Cray
          target_link_options( ectrans_lam_gpu_${prec} INTERFACE
                $<$<LINK_LANGUAGE:Fortran>:SHELL:${OpenMP_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:C,Cray>:SHELL:${OpenMP_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:CXX,Cray>:SHELL:${OpenMP_Fortran_FLAGS}> )
        endif()

        if( HAVE_ACC AND CMAKE_Fortran_COMPILER_ID MATCHES NVHPC )
          # Propagate flags as link options for downstream targets. Only required for NVHPC
          target_link_options( ectrans_lam_gpu_${prec} INTERFACE
                $<$<LINK_LANGUAGE:Fortran>:SHELL:${OpenACC_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:C,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS}>
                $<$<LINK_LANG_AND_ID:CXX,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS}> )
        endif()

        # This interface library is for backward compatibility, and provides the older includes
        ecbuild_add_library( TARGET etrans_gpu_${prec} TYPE INTERFACE )
        target_include_directories( etrans_gpu_${prec} INTERFACE $<BUILD_INTERFACE:${BUILD_INTERFACE_INCLUDE_DIR}/etrans_gpu_${prec}> )
        target_include_directories( etrans_gpu_${prec} INTERFACE $<INSTALL_INTERFACE:include/ectrans/etrans_gpu_${prec}> )
		target_link_libraries( etrans_gpu_${prec} INTERFACE ectrans_lam_gpu_${prec})


        # ## Install trans_gpu_${prec} interface
        # file( GLOB trans_interface ${PROJECT_SOURCE_DIR}/src/trans/include/ectrans/* )
        # install(
        #   FILES        ${trans_interface}
        #   DESTINATION  include/ectrans/trans_gpu_${prec}
        # )
  endif()
endforeach()

## Install etrans interface
install( DIRECTORY ${BUILD_INTERFACE_INCLUDE_DIR}/ DESTINATION include/ectrans )
