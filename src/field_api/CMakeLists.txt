#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

function(generate_file)
  set (options)
  set (oneValueArgs INPUT OUTPUT BACKEND)
  set (multiValueArgs)
  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(output  ${_PAR_OUTPUT})
  set(input   ${_PAR_INPUT})
  set(backend ${_PAR_BACKEND})
  set(sed_rules ${PROJECT_SOURCE_DIR}/src/field_api/sedrenames.txt)

  set( JPRB_dp JPRD )
  set( jprb_dp jprd )
  set( JPRB_sp JPRM )
  set( jprb_sp jprm )
  set( JPRB_gpu_dp JPRD )
  set( jprb_gpu_dp jprd )
  set( JPRB_gpu_sp JPRM )
  set( jprb_gpu_sp jprm )

  add_custom_command(
    OUTPUT  ${output}
    COMMAND  cat ${sed_rules} |
             sed -e "s/VARIANTDESIGNATOR/${backend}/g" |
             sed -e "s/TYPEDESIGNATOR_UPPER/${JPRB_${backend}}/g" |
             sed -e "s/TYPEDESIGNATOR_LOWER/${jprb_${backend}}/g" |
             sed -rf - ${input} > ${output}
    DEPENDS  ${input} ${sed_rules}
    COMMENT "Generating ${output}"
    VERBATIM
  )
endfunction(generate_file)

function(generate_backend_sources)
  set (options)
  set (oneValueArgs BACKEND DESTINATION OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})

  set(files)
  list( APPEND files
    field_api_ectrans_mod.F90
    trans/dir_trans_field_api.F90
    trans/inv_trans_field_api.F90
    field_api_basic_type_mod.F90
  )

  set(outfiles)
  foreach(file_i ${files})
    get_filename_component(outfile_name    ${file_i} NAME)
    get_filename_component(outfile_name_we ${file_i} NAME_WE)
    get_filename_component(outfile_ext     ${file_i} EXT)
    get_filename_component(outfile_dir     ${file_i} DIRECTORY)
    set(outfile "${destination}/${file_i}")

    generate_file(BACKEND ${backend} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile})
    list(APPEND outfiles ${outfile})
  endforeach(file_i)
  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)

set( BUILD_INTERFACE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include/ectrans )

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
 
      generate_backend_sources(BACKEND ${platform_tag}${prec}
                               OUTPUT ectrans_field_api_${platform_tag}${prec}_src
                               DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_field_api_${platform_tag}${prec})

      ecbuild_add_library(
          TARGET           ectrans_field_api_${platform_tag}${prec}
          LINKER_LANGUAGE  Fortran
          SOURCES          ${ectrans_field_api_${platform_tag}${prec}_src}
          PUBLIC_LIBS      ectrans_${platform_tag}${prec} field_api_${prec} fiat parkind_${prec}
          PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                           $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/field_api/include>
                           )

      ectrans_target_fortran_module_directory(
          TARGET            ectrans_field_api_${platform_tag}${prec}
          MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans
          INSTALL_DIRECTORY module/ectrans
        )

      if (${platform} MATCHES "gpu")
          if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
             set_source_files_properties(  ${ectrans_field_api_${platform_tag}${prec}_src} PROPERTIES COMPILE_OPTIONS "-acc" )
             #     set_source_files_properties(  ${ectrans_field_api_${platform_tag}${prec}_src} PROPERTIES COMPILE_OPTIONS "-O1" )
          endif()
      endif()
     endif()

  endforeach()
endforeach()



