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

foreach( program ectrans-benchmark ectrans-benchmark-ifs )

  if ( HAVE_CPU )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-${prec}
                                SOURCES ${program}.F90
                                LINKER_LANGUAGE Fortran
                                LIBS
                                  fiat
                                  parkind_${prec}
                                  trans_${prec}
        )
      endif()
    endforeach( prec)
  endif()

  if( HAVE_GPU )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        foreach( gpumethod acc omp )
          if( HAVE_${gpumethod} )
            ecbuild_add_executable( TARGET  ${program}-gpu-${prec}-${gpumethod}
                                    SOURCES ${program}.F90
                                    LINKER_LANGUAGE Fortran
                                    LIBS
                                      fiat
                                      parkind_${prec}
                                      trans_gpu_${prec}_${gpumethod}
            )

            if( gpumethod STREQUAL acc )
              target_link_libraries( ${program}-gpu-${prec}-${gpumethod} OpenACC::OpenACC_Fortran )
              target_compile_definitions( ${program}-gpu-${prec}-${gpumethod} PRIVATE ACCGPU )
            endif()

            if( gpumethod STREQUAL omp )
              target_compile_definitions( ${program}-gpu-${prec}-${gpumethod} PRIVATE OMPGPU )
            endif()

            target_link_options   ( ${program}-gpu-${prec}-${gpumethod} INTERFACE "$<$<LINK_LANG_AND_ID:Fortran,NVHPC>:${${gpumethod}_link_flags}>" )
            target_compile_options( ${program}-gpu-${prec}-${gpumethod} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${${gpumethod}_flags}>)

          endif()
        endforeach( gpumethod )
      endif()
    endforeach( prec )
  endif()
endforeach( program )

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
