# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

#foreach( program ectrans-benchmark ectrans-benchmark-ifs )
foreach( program ectrans-benchmark )

  list( APPEND util_src
    util/ectrans_memory.F90
    util/ectrans_memory.c
  )

  if ( HAVE_CPU )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-cpu-${prec}
                                SOURCES ${program}.F90 ${util_src}
                                LINKER_LANGUAGE Fortran
                                LIBS
                                  fiat
                                  parkind_${prec}
                                  trans_${prec}
                                DEFINITIONS
                                  VERSION="cpu"
        )
        ectrans_target_fortran_module_directory(TARGET ${program}-cpu-${prec} MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${program}-cpu-${prec} )
      endif()
    endforeach( prec)
  endif()

  if( HAVE_GPU )
    list( APPEND util_PRIVATE_DEFINITIONS
      $<${HAVE_CUDA}:CUDA>
      $<${HAVE_HIP}:HIP>
    )
    list( APPEND util_PRIVATE_LIBRARIES
      $<${HAVE_CUDA}:CUDA::cudart>
      $<${HAVE_HIP}:hip::host>
    )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-gpu-${prec}
                                SOURCES ${program}.F90 ${util_src}
                                LINKER_LANGUAGE Fortran
                                LIBS
                                  fiat
                                  parkind_${prec}
                                  trans_gpu_${prec}
                                  $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                                  $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                                  ${util_PRIVATE_LIBRARIES}
                                DEFINITIONS
                                  VERSION="gpu"
                                  ${util_PRIVATE_DEFINITIONS}
        )
        ectrans_target_fortran_module_directory(TARGET ${program}-gpu-${prec} MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${program}-gpu-${prec} )
      endif()
    endforeach( prec )
  endif( HAVE_GPU )
endforeach( program )

if( HAVE_ETRANS )
  foreach( program ectrans-lam-benchmark )

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
                    etrans_${prec}
                  DEFINITIONS
                    VERSION="cpu"
      )
      ectrans_target_fortran_module_directory(TARGET ${program}-cpu-${prec} MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${program}-cpu-${prec} )
      endif()
    endforeach( prec)
    endif()

    if( HAVE_GPU AND TRUE )
    list( APPEND util_PRIVATE_DEFINITIONS
      $<${HAVE_CUDA}:CUDA>
      $<${HAVE_HIP}:HIP>
    )
    list( APPEND util_PRIVATE_LIBRARIES
      $<${HAVE_CUDA}:CUDA::cudart>
      $<${HAVE_HIP}:hip::host>
    )
    foreach( prec dp sp )
      if( HAVE_${prec} )
      ecbuild_add_executable( TARGET  ${program}-gpu-${prec}
                  SOURCES ${program}.F90
                  LINKER_LANGUAGE Fortran
                  LIBS
                    fiat
                    parkind_${prec}
                    trans_gpu_${prec}
                    etrans_gpu_${prec}
                    $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                    $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                    ${util_PRIVATE_LIBRARIES}
                  DEFINITIONS
                    VERSION="gpu"
                    ${util_PRIVATE_DEFINITIONS}
      )
      ectrans_target_fortran_module_directory(TARGET ${program}-gpu-${prec} MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${program}-gpu-${prec} )
      endif()
    endforeach( prec )
    endif( )
  endforeach( program )
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
