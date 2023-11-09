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

ecbuild_list_add_pattern( LIST driver_src
                          GLOB
			    ectrans-benchmark-driver_mod.F90
			    ectrans-benchmark-data_mod.F90
                          QUIET
                        )
foreach( prec sp dp )
  if( HAVE_${prec} )
     if( ${prec} MATCHES "sp" )
       set(precno 1)
     elseif (${prec} MATCHES "dp")
       set(precno 2)
     endif()
    ecbuild_add_library(
      TARGET           ectrans_interfaces_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          wrappers/dir_trans.F90
                       wrappers/inv_trans.F90 
                       wrappers/setup_trans0.F90
                       wrappers/setup_trans.F90
                       wrappers/specnorm.F90
                       wrappers/trans_inq.F90
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
		       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PRIVATE_DEFINITIONS SYMBOLSUFFIX 
                       PRECOPT=${precno}
      PUBLIC_LIBS      fiat parkind_${prec} trans_dp trans_sp
    )

    ecbuild_add_executable(TARGET  ectrans-benchmark-${prec}
                           SOURCES ectrans-benchmark.F90
                           LIBS
                             fiat
                             parkind_${prec}
			       trans_${prec}
			       ectrans_interfaces_${prec}

                          )
  endif()
  if( HAVE_${prec} )
     if( ${prec} MATCHES "sp" )
       set(precno 1)
     elseif (${prec} MATCHES "dp")
       set(precno 2)
     endif()

  ecbuild_add_library(
      TARGET           driver_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${driver_src}
      PUBLIC_LIBS      fiat 
      PRIVATE_LIBS     parkind_${prec} trans_${prec}
   PRIVATE_DEFINITIONS SYMBOLSUFFIX=${prec}
                            PRECOPT=${precno}

    )

    ectrans_target_fortran_module_directory(
      TARGET            driver_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/driver_${prec}
      INSTALL_DIRECTORY module/driver_${prec}
    )

#   ecbuild_add_executable(TARGET  ectrans-benchmark-uni-${prec}
#                          SOURCES 
#       		   ectrans-benchmark-program.F90 
#                          LIBS
#                            fiat
#                            parkind_${prec}
#                            trans_${prec}
#                            driver_${prec}
#                            DEFINITIONS "SYMBOLSUFFIX=${prec}"
#                         )
 endif()
endforeach()
    ecbuild_add_executable(TARGET  ectrans-benchmark-uni
                           SOURCES 
			   ectrans-benchmark-program.F90 
                           LIBS
                             fiat
                             parkind_sp
                             trans_sp
                             driver_sp
                             parkind_dp
                             trans_dp
                             driver_dp   
                          )

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

