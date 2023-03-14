# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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
                          QUIET
                        )

## for reduced memory option, replace source files
if( HAVE_GPU_REDUCED_MEMORY )
  ecbuild_list_add_pattern( LIST reducedmem_files
                            GLOB internal_reducedmem/* QUIET )
  foreach( src_file ${reducedmem_files} )
    get_filename_component( base_name ${src_file} NAME )
    list(REMOVE_ITEM trans_src "internal/${base_name}" )
    list(APPEND      trans_src ${src_file} )
  endforeach()
endif()

ecbuild_list_exclude_pattern( LIST trans_src            REGEX dilatation_mod.F90 )

set_source_files_properties( internal/ftinv_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
set_source_files_properties( internal/ftdir_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )

ecbuild_list_add_pattern( LIST trans_hip_src
                          GLOB
                                sharedmem/*.hip.cpp
                                algor/*.hip.cpp
                                internal/*.hip.cpp
                                external/*.hip.cpp
                          QUIET
                        )
                ecbuild_warn_var(trans_hip_src)
                set_source_files_properties( ${trans_hip_src} PROPERTIES COMPILE_FLAGS "--offload-arch=gfx90a -x hip" )

foreach( prec sp )
  if( HAVE_${prec} )
    foreach( gpumethod acc)
      if( HAVE_${gpumethod})

        # We build an object library first. And then use these objects to create a shared and static library.

        ectrans_add_library(
          TARGET           trans_gpu_object_${prec}_${gpumethod}
          TYPE             OBJECT
          SOURCES          ${trans_src}
          PUBLIC_INCLUDES  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
                          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/ectrans>
                          $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/algor/interface>
                          $<INSTALL_INTERFACE:include/ectrans>
                          $<INSTALL_INTERFACE:include>
          PUBLIC_LIBS      parkind_${prec}
                          fiat
          PRIVATE_LIBS    hip::hipfft
                          roc::hipblas
			                    roc::rocblas
			                    roc::rocfft
                          MPI::MPI_Fortran
                          mpifort
                          omptarget
                          #MPI::MPI_C
			  MPI::MPI_CXX
                          #nvhpcwrapnvtx
                          #${LAPACK_LIBRARIES}
        )

        target_link_options   ( trans_gpu_object_${prec}_${gpumethod} INTERFACE "$<$<LINK_LANG_AND_ID:Fortran,NVHPC>:${${gpumethod}_link_flags}>" )
        target_compile_options( trans_gpu_object_${prec}_${gpumethod} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${${gpumethod}_flags}>)

        ectrans_target_fortran_module_directory(
          TARGET            trans_gpu_object_${prec}_${gpumethod}
          MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_gpu_${prec}_${gpumethod}
          INSTALL_DIRECTORY module/trans_gpu_${prec}_${gpumethod}
        )

        if( HAVE_OMP )
          target_link_libraries( trans_gpu_object_${prec}_${gpumethod} PRIVATE OpenMP::OpenMP_Fortran )
        endif()

        if( prec STREQUAL sp )
          target_compile_definitions( trans_gpu_object_${prec}_${gpumethod} PRIVATE TRANS_SINGLE PARKINDTRANS_SINGLE )
        endif()

        if( gpumethod STREQUAL acc )
          target_compile_definitions( trans_gpu_object_${prec}_${gpumethod} PRIVATE ACCGPU )
          target_link_libraries( trans_gpu_object_${prec}_${gpumethod} PRIVATE OpenACC::OpenACC_Fortran )
        endif()

        if( gpumethod STREQUAL omp )
          target_compile_definitions( trans_gpu_object_${prec}_${gpumethod} PRIVATE OMPGPU )
        endif()

        if( HAVE_GPU_REDUCED_MEMORY )
          target_compile_definitions( trans_gpu_object_${prec}_${gpumethod} PRIVATE REDUCED_MEM )
        endif()

        if( HAVE_GPU_AWARE_MPI )
          target_compile_definitions( trans_gpu_object_${prec}_${gpumethod} PRIVATE USE_CUDA_AWARE_MPI_FT )
        endif()

        ectrans_add_library(
          TARGET           trans_gpu_shared_${prec}_${gpumethod}
          OUTPUT_NAME      trans_gpu_${prec}_${gpumethod}
          TYPE             SHARED
          LINKER_LANGUAGE  Fortran
          PUBLIC_LIBS      trans_gpu_object_${prec}_${gpumethod}
        )

        ectrans_add_library(
          TARGET           trans_gpu_${prec}_${gpumethod}
          TYPE             STATIC
          LINKER_LANGUAGE  Fortran
          PUBLIC_LIBS      trans_gpu_object_${prec}_${gpumethod}
        )
      endif()
    endforeach()
  endif()
endforeach()
