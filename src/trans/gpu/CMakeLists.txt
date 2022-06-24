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
      FLAGS   "-fast -O3 -acc -gpu=pinned,deepcopy,fastmath,nordc")

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


## for reduced memory option, remove duplicate source files
if( HAVE_REDUCED_MEMORY )
  ecbuild_list_add_pattern( LIST reducedmem_files
                            GLOB internal_reducedmem/* QUIET )
  ecbuild_list_add_pattern( LIST trans_src
	                    GLOB internal_reducedmem/* QUIET )
  foreach( src_file ${reducedmem_files} )
    get_filename_component( base_name ${src_file} NAME )
    list(REMOVE_ITEM trans_src "internal/${base_name}" )
  endforeach()
endif()

ecbuild_list_exclude_pattern( LIST trans_src REGEX dilatation_mod.F90 )
ecbuild_list_exclude_pattern( LIST trans_src_reducedmem REGEX dilatation_mod.F90 )

set( Debug_Compile_Flags            -g -Mlarge_arrays -traceback )
set( Release_Compile_Flags       -fast -Mlarge_arrays )
set( RelDeb_Compile_Flags  -gopt -fast -Mlarge_arrays )

foreach( prec sp dp )
  if( HAVE_${prec} )

    # We build an object library first. And then use these objects to create a shared and static library.

    ectrans_add_library(
      TARGET           trans_gpu_object_${prec}
      TYPE             OBJECT
      SOURCES          ${trans_src}
                       algor/external/fourier/destroy_plan_fftc.cu
                       algor/external/fourier/create_plan_fftc.cu
                       algor/external/fourier/storage_fftc.cu
                       algor/external/fourier/execute_plan_fftc.cu
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/ectrans>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/algor/interface>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PRIVATE_INCLUDES ${MPI_Fortran_INCLUDE_PATH}
      PUBLIC_LIBS      parkind_${prec}
                       fiat
      PRIVATE_LIBS     ${LAPACK_LIBRARIES} OpenACC::OpenACC_Fortran
    )
    ectrans_target_fortran_module_directory(
      TARGET            trans_gpu_object_${prec}
      MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_gpu_${prec}
      INSTALL_DIRECTORY module/trans_gpu_${prec}
    )

    if( HAVE_OMP )
      target_link_libraries( trans_gpu_object_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    target_compile_options( trans_gpu_object_${prec} PRIVATE
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:${Debug_Compile_Flags}>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:${Release_Compile_Flags}>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RelWithDebInfo>>:${RelDeb_Compile_Flags}>" )

    if( prec STREQUAL sp )
      target_compile_definitions( trans_gpu_object_${prec} PRIVATE TRANS_SINGLE PARKINDTRANS_SINGLE )
    endif()
    if( HAVE_REDUCED_MEMORY )
      target_compile_definitions( trans_gpu_object_${prec} PRIVATE REDUCED_MEM )
    endif()

    set( HAVE_CUDA_AWARE_MPI 1 )
    if( HAVE_CUDA_AWARE_MPI )
      target_compile_definitions( trans_gpu_object_${prec} PRIVATE USE_CUDA_AWARE_MPI_FT )
    endif()


    ectrans_add_library(
      TARGET           trans_gpu_shared_${prec}
      OUTPUT_NAME      trans_gpu_${prec}
      TYPE             SHARED
      LINKER_LANGUAGE  Fortran
      PUBLIC_LIBS      trans_gpu_object_${prec}
    )

    ectrans_add_library(
      TARGET           trans_gpu_${prec}
      TYPE             STATIC
      LINKER_LANGUAGE  Fortran
      PUBLIC_LIBS      trans_gpu_object_${prec}
    )
  endif()
endforeach()

## precision-independent GPU library with CUDA kernels
ecbuild_add_library( TARGET gpu_cuda
  TYPE STATIC
  SOURCES
            algor/module/cublasSgemmBatched.cu
            algor/module/cublasDgemmBatched.cu
            algor/module/IPC_Alltoall.cu
  PRIVATE_LIBS MPI::MPI_C
)
