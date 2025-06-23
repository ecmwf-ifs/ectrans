# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Sources which are precision independent can go into a common library
list( APPEND ectrans_lam_common_src
  internal/tpmald_distr.F90
  internal/tpmald_dim.F90
  internal/esetup_dims_mod.F90
  internal/esetup_geom_mod.F90
  internal/suemplat_mod.F90
  internal/suemplatb_mod.F90
  internal/ellips.F90
)
list( APPEND ectrans_lam_common_includes
)

ecbuild_add_library(
  TARGET           ectrans_lam_common
  LINKER_LANGUAGE  Fortran
  SOURCES          ${ectrans_lam_common_src}
  PUBLIC_LIBS      fiat ectrans_common
  PRIVATE_LIBS     ${LAPACK_LIBRARIES}
  PUBLIC_INCLUDES  $<INSTALL_INTERFACE:include/ectrans>
                   $<INSTALL_INTERFACE:include>
                   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
)
ectrans_target_fortran_module_directory(
  TARGET            ectrans_lam_common
  MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans
  INSTALL_DIRECTORY module/ectrans
)

if( HAVE_OMP )
  ecbuild_debug("target_link_libraries( ectrans_common PRIVATE OpenMP::OpenMP_Fortran )")
  target_link_libraries( ectrans_lam_common PRIVATE OpenMP::OpenMP_Fortran )
endif()

set( ectrans_lam_common_includes ${ectrans_lam_common_includes} PARENT_SCOPE )
