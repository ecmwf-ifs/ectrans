# (C) Copyright 2026- ECMWF.
# (C) Copyright 2026- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


function(generate_backend_sources)
  set (options)
  set (oneValueArgs BACKEND DESTINATION PLATFORM OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})
  set(platform ${_PAR_PLATFORM})

  file(MAKE_DIRECTORY ${destination}/field_view)

  list( APPEND files
        ectrans_field_view_internal_util_mod.F90
        dir_trans_field_view.F90
        inv_trans_field_view.F90
        ${platform}/internal/dir_trans_view_ctl_mod.F90
        ${platform}/internal/inv_trans_view_ctl_mod.F90
        ${platform}/internal/ltdir_view_mod.F90
        ${platform}/internal/ltinv_view_mod.F90
        ${platform}/internal/trltog_view_mod.F90
        ${platform}/internal/trgtol_view_mod.F90
        ${platform}/internal/updsp_view_mod.F90
        ${platform}/internal/updspb_view_mod.F90
        ${platform}/internal/prfi1b_view_mod.F90)

  if (${platform} STREQUAL "cpu")
      list( APPEND files
        ${platform}/internal/trgl_view_mod.F90
        ${platform}/internal/fsc_view_mod.F90        
        ${platform}/internal/ftdir_view_ctl_mod.F90
        ${platform}/internal/ftinv_view_ctl_mod.F90
        ${platform}/internal/ltdir_view_ctl_mod.F90
        ${platform}/internal/ltinv_view_ctl_mod.F90)
  endif()

  set(outfiles)

  foreach(file_i ${files})
    set(outfile "${destination}/${file_i}")
    ecbuild_debug("Generate ${outfile}")
    generate_file(BACKEND ${backend} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile})
    list(APPEND outfiles ${outfile})
  endforeach(file_i)

  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)


set (platforms)
if (HAVE_CPU)
    list(APPEND platforms cpu)
endif()

# We can't reliably pass GPU buffers directly to MPL routines yet, so for now we should still use raw MPI calls
if( HAVE_GPU_AWARE_MPI )
  set( USE_RAW_MPI 1 )
else()
  set( USE_RAW_MPI 0)
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

      # Make sure the TYPE of ectrans_field_view respects HAVE_GPU_STATIC
      set( LIBRARY_TYPE SHARED )
      if( ${platform} MATCHES "gpu" AND HAVE_GPU_STATIC )
        set( LIBRARY_TYPE STATIC )
      endif()

      generate_backend_sources( BACKEND ${platform_tag}${prec}
        PLATFORM ${platform}
        OUTPUT ectrans_field_view_${platform_tag}${prec}_src
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_${platform_tag}${prec})
      ecbuild_add_library(
          TARGET           ectrans_field_view_${platform_tag}${prec}
          TYPE             ${LIBRARY_TYPE}
          LINKER_LANGUAGE  Fortran
          SOURCES          ${ectrans_field_view_${platform_tag}${prec}_src}
                           ectrans_field_view_mod.F90
          PUBLIC_LIBS      trans_${platform_tag}${prec} fiat parkind_${prec}
          PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/field_view/include/ectrans>
                           $<INSTALL_INTERFACE:include/ectrans> 
          PRIVATE_LIBS         ${ECTRANS_GPU_HIP_LIBRARIES}
                               $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                               $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                               $<${USE_RAW_MPI}:MPI::MPI_Fortran>
          PRIVATE_DEFINITIONS  ${GPU_RUNTIME}GPU ${GPU_OFFLOAD}GPU
                               $<${HAVE_CUTLASS}:USE_CUTLASS>
                               $<${HAVE_CUTLASS_3XTF32}:USE_CUTLASS_3XTF32>
                               $<${HAVE_GPU_GRAPHS_GEMM}:USE_GRAPHS_GEMM>
                               $<${HAVE_GPU_GRAPHS_FFT}:USE_GRAPHS_FFT>
                               $<${HAVE_ROCFFT}:USE_ROCFFT>
                               $<${HAVE_GPU_AWARE_MPI}:USE_GPU_AWARE_MPI>
                               $<${USE_RAW_MPI}:USE_RAW_MPI>
                             )
      ecbuild_target_fortran_module_directory(
          TARGET            ectrans_field_view_${platform_tag}${prec}
          MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans_field_view_${platform_tag}${prec}
          INSTALL_DIRECTORY module/ectrans_field_view_${platform_tag}${prec} )
          if( HAVE_OMP )
            ecbuild_debug("target_link_libraries( ectrans_field_view_${platform_tag}${prec} PRIVATE OpenMP::OpenMP_Fortran )")
            target_link_libraries( ectrans_field_view_${platform_tag}${prec} PRIVATE OpenMP::OpenMP_Fortran )
          endif()
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
