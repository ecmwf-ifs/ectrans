# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( NOT ectrans_VERSION_PATCH )
  set( ectrans_VERSION_PATCH 0 )
endif()

configure_file( version.c.in  version.c )

ecbuild_add_library( TARGET transi_dp
  SOURCES             transi_module.F90
                      transi.h
                      transi.c
                      version.h
                      ${CMAKE_CURRENT_BINARY_DIR}/version.c
  HEADER_DESTINATION  include/ectrans
  PUBLIC_INCLUDES     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                      $<INSTALL_INTERFACE:include>
  PRIVATE_LIBS        trans_dp
  PRIVATE_DEFINITIONS ECTRANS_HAVE_MPI=${ectrans_HAVE_MPI}
)
ectrans_target_fortran_module_directory( TARGET transi_dp
  MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module/trans_dp
)

if( HAVE_GPU )
  ecbuild_add_library( TARGET transi_gpu_dp
    SOURCES             transi_module.F90
                        transi.h
                        transi.c
                        version.h
                        ${CMAKE_CURRENT_BINARY_DIR}/version.c
    HEADER_DESTINATION  include/ectrans
    PUBLIC_INCLUDES     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                        $<INSTALL_INTERFACE:include>
    PRIVATE_LIBS        trans_gpu_dp
    PRIVATE_DEFINITIONS ECTRANS_HAVE_MPI=${ectrans_HAVE_MPI} ECTRANS_GPU_VERSION
  )
  if( HAVE_ACC AND CMAKE_Fortran_COMPILER_ID MATCHES NVHPC )
    # Propagate flags as link options for downstream targets. Only required for NVHPC
    target_link_options( transi_gpu_dp INTERFACE
      $<$<LINK_LANGUAGE:Fortran>:SHELL:${OpenACC_Fortran_FLAGS}>
      $<$<LINK_LANG_AND_ID:C,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS}>
      $<$<LINK_LANG_AND_ID:CXX,NVHPC>:SHELL:${OpenACC_Fortran_FLAGS}> )
  endif()

  ectrans_target_fortran_module_directory( TARGET transi_gpu_dp
    MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module/trans_gpu_dp
  )
endif()

set( transi_includes transi.h version.h )
install(
  FILES        ${transi_includes}
  DESTINATION  include/ectrans
)
