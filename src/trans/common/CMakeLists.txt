# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Preprocess module file containing version information
configure_file( internal/ectrans_version_mod.F90.in internal/ectrans_version_mod.F90 )

## Sources which are precision independent can go into a common library
list( APPEND ectrans_common_src
  sharedmem/sharedmem_mod.F90
  sharedmem/sharedmem.c
  internal/ectrans_blas_mod.F90
  internal/abort_trans_mod.F90
  internal/cpledn_mod.F90
  internal/field_split_mod.F90
  internal/gawl_mod.F90
  internal/sugaw_mod.F90
  internal/supol_mod.F90
  internal/supolf_mod.F90
  internal/tpm_constants.F90
  internal/tpm_ctl.F90
  internal/tpm_dim.F90
  internal/tpm_gen.F90
  internal/tpm_geometry.F90
  internal/tpm_pol.F90
  internal/tpm_distr.F90
  internal/pe2set_mod.F90
  internal/set2pe_mod.F90
  internal/eq_regions_mod.F90
  internal/shuffle_mod.F90
  internal/sump_trans0_mod.F90
  internal/sustaonl_mod.F90
  internal/sumplat_mod.F90
  internal/sumplatb_mod.F90
  internal/sumplatbeq_mod.F90
  internal/sumplatf_mod.F90
  internal/mysendset_mod.F90
  internal/myrecvset_mod.F90
  internal/suwavedi_mod.F90
  internal/sump_trans_preleg_mod.F90
  external/get_current.F90
  external/setup_trans0.F90
  external/ini_spec_dist.F90
  ${CMAKE_CURRENT_BINARY_DIR}/internal/ectrans_version_mod.F90
)
list( APPEND ectrans_common_includes
  ectrans/get_current.h
  ectrans/setup_trans0.h
  ectrans/ini_spec_dist.h
)

ecbuild_add_library(
  TARGET           ectrans_common
  LINKER_LANGUAGE  Fortran
  SOURCES          ${ectrans_common_src}
  PUBLIC_LIBS      fiat
  PRIVATE_LIBS     ${LAPACK_LIBRARIES}
  PUBLIC_INCLUDES  $<INSTALL_INTERFACE:include/ectrans>
                   $<INSTALL_INTERFACE:include>
                   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
)
ectrans_target_fortran_module_directory(
  TARGET            ectrans_common
  MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans
  INSTALL_DIRECTORY module/ectrans
)

if( HAVE_OMP )
  ecbuild_debug("target_link_libraries( ectrans_common PRIVATE OpenMP::OpenMP_Fortran )")
  target_link_libraries( ectrans_common PRIVATE OpenMP::OpenMP_Fortran )
endif()

set( ectrans_common_includes ${ectrans_common_includes} PARENT_SCOPE )
