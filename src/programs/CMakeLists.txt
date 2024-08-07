# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

#foreach( program ectrans-benchmark ectrans-benchmark-ifs )
foreach( program ectrans-benchmark )

  if ( HAVE_CPU )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-cpu-${prec}
                                SOURCES ${program}.F90
                                LINKER_LANGUAGE Fortran
                                LIBS
                                  fiat
                                  parkind_${prec}
                                  trans_${prec}
        )
        target_compile_definitions(${program}-cpu-${prec} PRIVATE VERSION="cpu")
      endif()
    endforeach( prec)
  endif()

  if( HAVE_GPU )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-gpu-${prec}
                                SOURCES ${program}.F90
                                LINKER_LANGUAGE Fortran
                                LIBS
                                  fiat
                                  parkind_${prec}
                                  trans_gpu_${prec}
                                  OpenACC::OpenACC_Fortran
       )
       target_compile_definitions(${program}-gpu-${prec} PRIVATE VERSION="gpu")
       endif()
    endforeach( prec )
  endif( HAVE_GPU )
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
