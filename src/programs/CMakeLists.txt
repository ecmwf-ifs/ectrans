# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( HAVE_TOOLS AND TARGET eccodes_f90 )

  if( HAVE_DOUBLE_PRECISION )
    set( trans trans_dp )
  else()
    set( trans trans_sp )
  endif()

  ecbuild_warn_var( ECTRANS_TOOLS_RTABLE_PATH )

  file( GLOB ectrans_programs *.F90 )
  foreach( _program IN ITEMS ${ectrans_programs} )
    get_filename_component( _program ${_program} NAME_WE )

    ecbuild_add_executable(TARGET ${_program}
      SOURCES ${_program}.F90
      LIBS ${trans} eccodes_f90
      LINKER_LANGUAGE Fortran
      DEFINITIONS ECTRANS_TOOLS_RTABLE_PATH="${ECTRANS_TOOLS_RTABLE_PATH}" )

  endforeach()


endif()


set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

if( HAVE_GPU )
  foreach( prec sp dp )
    if( HAVE_${prec} )
      ecbuild_add_executable(TARGET  driver-spectrans-${prec}
                             SOURCES driver-spectraltransform.F90
                             INCLUDES
                                  ${MPI_Fortran_INCLUDE_PATH}
                                  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../trans/gpu/include/ectrans>
                             LIBS 
                                  fiat parkind_${prec}
                                  eccodes_f90 eccodes_memfs
                                  ${MPI_Fortran_LIBRARIES}
                                  trans_gpu_static_${prec}
                                  gpu
                                  OpenACC::OpenACC_Fortran 
                                  ${LAPACK_LIBRARIES}
                                  nvhpcwrapnvtx
                          )
	ecbuild_add_executable(TARGET  driver-spectrans-CA-${prec}
                             SOURCES driver-spectraltransform.F90
                             INCLUDES
                                  ${MPI_Fortran_INCLUDE_PATH}
                                  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../trans/gpu/include/ectrans>
                             LIBS 
                                  fiat parkind_${prec}
                                  eccodes_f90 eccodes_memfs
                                  ${MPI_Fortran_LIBRARIES}
                                  trans_gpu_static_CA_${prec}
                                  gpu
                                  OpenACC::OpenACC_Fortran 
                                  ${LAPACK_LIBRARIES}
                                  nvhpcwrapnvtx
                          )
		  #trans_gpu_static_${prec}
		  #gpu
		  #${CMAKE_BINARY_DIR}/lib/libtrans_gpu_static_${prec}.a
		  #${CMAKE_BINARY_DIR}/lib/libgpu.a
      #target_link_libraries( driver-spectrans PRIVATE OpenACC::OpenACC_Fortran )
      set_property( TARGET driver-spectrans-${prec} PROPERTY CUDA_ARCHITECTURES 70 )
      set_property( TARGET driver-spectrans-CA-${prec} PROPERTY CUDA_ARCHITECTURES 70 )
      target_compile_options( driver-spectrans-${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-g -acc -Minfo=acc -gpu=cc70,lineinfo,deepcopy,fastmath,nordc -cudalib=cufft,cublas -fpic>  )
      target_compile_options( driver-spectrans-CA-${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-g -acc -Minfo=acc -gpu=cc70,lineinfo,deepcopy,fastmath,nordc -cudalib=cufft,cublas -fpic>  )
      set_target_properties(driver-spectrans-${prec} PROPERTIES LINK_FLAGS "-acc -cudalib=cufft,cublas -fpic")
      set_target_properties(driver-spectrans-CA-${prec} PROPERTIES LINK_FLAGS "-acc -cudalib=cufft,cublas -fpic")
      message("Building ${prec} GPU driver")
    endif()
  endforeach()
else()
  ecbuild_warn("HAVE_GPU false so not building GPU driver")
endif()

# ectrans information tool

get_property( langs GLOBAL PROPERTY ENABLED_LANGUAGES )

foreach( lang ${langs} )
  set( EC_${lang}_FLAGS "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}" )
endforeach()

configure_file( ectrans.in ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ectrans @ONLY )

file(COPY ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ectrans
  DESTINATION ${CMAKE_BINARY_DIR}/bin
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

install( FILES
  ${CMAKE_BINARY_DIR}/bin/ectrans
  DESTINATION
  ${INSTALL_BIN_DIR}
  PERMISSIONS
  OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

