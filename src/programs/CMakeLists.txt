# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( HAVE_TOOLS )

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
      LIBS ${trans}
      LINKER_LANGUAGE Fortran
      DEFINITIONS ECTRANS_TOOLS_RTABLE_PATH="${ECTRANS_TOOLS_RTABLE_PATH}" )

  endforeach()

endif()


set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

foreach( prec sp dp )
  if( HAVE_${prec} )
    ecbuild_add_executable(TARGET  ectrans-benchmark-${prec}
                           SOURCES ectrans-benchmark.F90
                           LIBS
                             fiat
                             parkind_${prec}
                             trans_${prec}
                          )
  endif()
endforeach()

if( HAVE_GPU )
  set( Debug_Compile_Flags   -g -O0 -acc -Minfo=acc -gpu=lineinfo,deepcopy,fastmath,nordc -cudalib=cufft,cublas -fpic -traceback )
  set( Release_Compile_Flags -fast -acc -Minfo=acc -gpu=lineinfo,deepcopy,fastmath,nordc -cudalib=cufft,cublas -fpic )
  set( RelDeb_Compile_Flags -gopt -fast -acc -Minfo=acc -gpu=lineinfo,deepcopy,fastmath,nordc -Mlarge_arrays )
  foreach( prec sp dp )
    if( HAVE_${prec} )
      ecbuild_add_executable(TARGET  ectrans-benchmark-gpu-${prec}
                             SOURCES ectrans-benchmark.F90
                             LIBS
                                  fiat
                                  parkind_${prec}
                                  trans_gpu_static_${prec}
                                  gpu
                                  OpenACC::OpenACC_Fortran
                                  nvhpcwrapnvtx
                          )
	ecbuild_add_executable(TARGET  ectrans-benchmark-gpu-CA-${prec}
                             SOURCES ectrans-benchmark.F90
                             LIBS
                                  fiat
                                  parkind_${prec}
                                  trans_gpu_static_CA_${prec}
                                  gpu
                                  OpenACC::OpenACC_Fortran
                                  ${LAPACK_LIBRARIES}
                                  nvhpcwrapnvtx
                          )
      target_compile_options( ectrans-benchmark-gpu-${prec} PRIVATE
	    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:${Debug_Compile_Flags}>"
	    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:${Release_Compile_Flags}>"
	    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RelWithDebInfo>>:${RelDeb_Compile_Flags}>" )
      target_compile_options( ectrans-benchmark-gpu-CA-${prec} PRIVATE
	    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:${Debug_Compile_Flags}>"
	    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:${Release_Compile_Flags}>"
	    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RelWithDebInfo>>:${RelDeb_Compile_Flags}>" )
      set_target_properties(ectrans-benchmark-gpu-${prec} PROPERTIES LINK_FLAGS "-acc -cudalib=cufft,cublas -fpic")
      set_target_properties(ectrans-benchmark-gpu-CA-${prec} PROPERTIES LINK_FLAGS "-acc -cudalib=cufft,cublas -fpic")
      message("Building ${prec} GPU benchmark")
    endif()
  endforeach()
else()
  ecbuild_warn("HAVE_GPU false so not building GPU benchmark")
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
