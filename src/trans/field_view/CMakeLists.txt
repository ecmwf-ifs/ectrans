# (C) Copyright 2026- ECMWF.
# (C) Copyright 2026- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set (platforms)
if (HAVE_CPU)
    list(APPEND platforms cpu)
endif()
if (HAVE_GPU)
    list(APPEND platforms gpu)
endif()
foreach(platform ${platforms})
  foreach( prec dp sp)
    if( HAVE_${prec})
      if (${platform} MATCHES "gpu")
        set (platform_tag "gpu_")
      else()
        set (platform_tag "")
      endif()

      ecbuild_add_library(
          TARGET           ectrans_field_view_${platform_tag}${prec}
          LINKER_LANGUAGE  Fortran
          SOURCES          ectrans_field_view_mod.F90
                           ectrans_field_view_internal_util_mod.F90
                           dir_trans_field_view.F90
                           inv_trans_field_view.F90
          PUBLIC_LIBS      trans_${platform_tag}${prec} fiat parkind_${prec}
          PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/field_view/include/ectrans>
                           $<INSTALL_INTERFACE:include/ectrans> )

      ecbuild_target_fortran_module_directory(
          TARGET            ectrans_field_view_${platform_tag}${prec}
          MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans_field_view_${platform_tag}${prec}
          INSTALL_DIRECTORY module/ectrans_field_view_${platform_tag}${prec} )
    endif()
  endforeach()
endforeach()

list( APPEND field_view_includes 
  include/ectrans/dir_trans_field_view.h
  include/ectrans/inv_trans_field_view.h
)
install(
  FILES        ${field_view_includes}
  DESTINATION  include/ectrans
)
