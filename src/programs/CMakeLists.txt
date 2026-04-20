# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

#foreach( program ectrans-benchmark ectrans-benchmark-ifs )
foreach( program ectrans-benchmark )

  ## TODO: Improve ectrans_checksum_mod.F90 to be precision independent (independent of parkind1 module).
  #        Since this module uses trans functions like gath_grid, gath_spec,
  #        this can only be made independent of precision when we have a dispatch library that can call
  #        the correct version of these functions for the precision of the input data.
  #        For now, we include this module in the ectrans-benchmark programs.
  list( APPEND util_src
    util/ectrans_checksum_mod.F90
  )

  ## NOTE: ectrans_field_api_helper is precision-dependent but only needed in ectrans-benchmark for now.
  #        It can only be made into a precision-independent module if field_api itself can
  #        be made independent of parkind1, not defining fields like FIELD_3RB (needing JPRB).
  if (HAVE_FIELD_API)
    list(APPEND util_src util/ectrans_field_api_helper.F90)
  endif()

  if ( HAVE_CPU )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-cpu-${prec}
                                SOURCES ${program}.F90 ${util_src}
                                LINKER_LANGUAGE Fortran
                                LIBS
                                  fiat
                                  parkind_${prec}
                                  ectrans_util_cpu
                                  trans_${prec}
                                  $<${HAVE_FIELD_API}:ectrans_field_api_${prec}>
                                DEFINITIONS
                                  VERSION="cpu"
                                  USE_FIELD_API=${HAVE_FIELD_API}
        )
        ecbuild_target_fortran_module_directory(TARGET ${program}-cpu-${prec} MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${program}-cpu-${prec} )
      endif()
    endforeach( prec)
  endif()

  if( HAVE_GPU )
    list( APPEND util_PRIVATE_DEFINITIONS
      $<${HAVE_CUDA}:CUDA>
      $<${HAVE_HIP}:HIP>
    )
    foreach( prec dp sp )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-gpu-${prec}
                                SOURCES ${program}.F90 ${util_src}
                                LINKER_LANGUAGE Fortran
                                LIBS
                                  fiat
                                  parkind_${prec}
                                  ectrans_util
                                  trans_gpu_${prec}
                                  $<${HAVE_FIELD_API}:ectrans_field_api_gpu_${prec}>
                                DEFINITIONS
                                  VERSION="gpu"
                                  USE_FIELD_API=${HAVE_FIELD_API}
        )
        ecbuild_target_fortran_module_directory(TARGET ${program}-gpu-${prec} MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${program}-gpu-${prec} )
      endif()
    endforeach( prec )
  endif( HAVE_GPU )
endforeach( program )

if( HAVE_ETRANS )
  set(program ectrans-lam-benchmark)
  foreach( prec sp dp )
    if( HAVE_${prec} )
      ecbuild_add_executable( TARGET  ${program}-cpu-${prec}
                              SOURCES ${program}.F90
                              LINKER_LANGUAGE Fortran
                              LIBS
                                fiat
                                parkind_${prec}
                                ectrans_util_cpu
                                trans_${prec}
                                etrans_${prec}
                                $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                            )
    endif()
  
    if( HAVE_ETRANS_GPU )
      if( HAVE_${prec} )
        ecbuild_add_executable( TARGET  ${program}-gpu-${prec}
                    SOURCES ${program}.F90
                    LINKER_LANGUAGE Fortran
                    LIBS
                      fiat
                      parkind_${prec}
                      ectrans_util
                      trans_gpu_${prec}
                      etrans_gpu_${prec}
                    DEFINITIONS
                      VERSION="gpu"
        )
        ecbuild_target_fortran_module_directory(TARGET ${program}-gpu-${prec} MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/${program}-gpu-${prec} )
      endif()
    endif( )
  endforeach( prec )
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
