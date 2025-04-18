# (C) Copyright 2020- ECMWF.
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
  set(sed_rules ${PROJECT_SOURCE_DIR}/src/trans/sedrenames.txt)

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


function(generate_backend_includes)
  set (options)
  set (oneValueArgs BACKEND TARGET DESTINATION INCLUDE_DIRECTORY)
  set (multiValueArgs)
  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(destination ${_PAR_DESTINATION} )
  set(backend     ${_PAR_BACKEND})

  file(MAKE_DIRECTORY ${destination})
  file(MAKE_DIRECTORY ${destination}/trans_${backend})

  ecbuild_list_add_pattern( LIST absolute_files GLOB ectrans/*.h SOURCE_DIR ${_PAR_INCLUDE_DIRECTORY} QUIET )
  set( files )
  foreach(file_i ${absolute_files})
    file(RELATIVE_PATH file_i ${_PAR_INCLUDE_DIRECTORY} ${file_i})
    list(APPEND files ${file_i})
  endforeach()
  set( outfiles )
  foreach(file_i ${files})
    get_filename_component(outfile_name    ${file_i} NAME)
    get_filename_component(outfile_name_we ${file_i} NAME_WE)
    get_filename_component(outfile_ext     ${file_i} EXT)
    get_filename_component(outfile_dir     ${file_i} DIRECTORY)
    if (${file_i} IN_LIST ectrans_common_includes)
      configure_file(${_PAR_INCLUDE_DIRECTORY}/${file_i} ${destination}/${outfile_name})
    else()
      set(outfile "${destination}/${outfile_name_we}_${backend}${outfile_ext}")
      ecbuild_debug("Generate ${outfile}")
      generate_file(BACKEND ${backend} INPUT ${_PAR_INCLUDE_DIRECTORY}/${file_i} OUTPUT ${outfile})
      list(APPEND outfiles ${outfile})
      string(TOUPPER ${outfile_name_we} OUTFILE_NAME_WE )
      ecbuild_debug("Generate ${destination}/trans_${backend}/${outfile_name}")
      file(WRITE  ${destination}/trans_${backend}/${outfile_name} "! Automatically generated interface header for backward compatibility of generic symbols !\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#if defined(${outfile_name_we})\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#undef ${outfile_name_we}\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#endif\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#if defined(${OUTFILE_NAME_WE})\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#undef ${OUTFILE_NAME_WE}\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#endif\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#include \"${outfile_name_we}_${backend}${outfile_ext}\"\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#define ${outfile_name_we} ${OUTFILE_NAME_WE}_${backend}\n")
      file(APPEND ${destination}/trans_${backend}/${outfile_name} "#define ${OUTFILE_NAME_WE} ${OUTFILE_NAME_WE}_${backend}\n")
    endif()
  endforeach(file_i)

  add_custom_target(${_PAR_TARGET}_generate DEPENDS ${outfiles})
  ecbuild_add_library(TARGET ${_PAR_TARGET} TYPE INTERFACE)
  add_dependencies(${_PAR_TARGET} ${_PAR_TARGET}_generate)
  target_include_directories(${_PAR_TARGET} INTERFACE $<BUILD_INTERFACE:${destination}>)
endfunction(generate_backend_includes)






add_subdirectory( common )

if( HAVE_CPU)
  add_subdirectory( cpu )
endif()
if( HAVE_GPU )
  add_subdirectory( gpu )
endif()
