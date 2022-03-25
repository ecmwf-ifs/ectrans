# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                algor/*
                                internal/*
                                external/*
                          QUIET
                        )

if( NOT HAVE_FFTW )
  ecbuild_list_exclude_pattern( LIST trans_src REGEX tpm_fftw.F90 )
endif()

set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )
set( FFTW_dp FFTW::fftw3  )
set( FFTW_sp FFTW::fftw3f )

foreach( prec sp dp )
  if( HAVE_${prec} )

    ecbuild_add_library( 
      TARGET           trans_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_src}
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ectrans>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PUBLIC_LIBS      parkind_${prec}
                       fiat
      PRIVATE_LIBS     ${LAPACK_LIBRARIES}
    )
    ectrans_target_fortran_module_directory(
      TARGET            trans_${prec}
      MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_${prec}
      INSTALL_DIRECTORY module/trans_${prec}
    )
    if( HAVE_FFTW )
      target_link_libraries( trans_${prec} PRIVATE ${FFTW_LIBRARIES} )
      target_include_directories( trans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
      target_compile_definitions( trans_${prec} PRIVATE WITH_FFTW )
    endif()
    if( HAVE_OMP )
      target_link_libraries( trans_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()
  endif()
endforeach()

## Install trans interface

file( GLOB trans_interface include/ectrans/* )
install(
  FILES        ${trans_interface}
  DESTINATION  include/ectrans
)
