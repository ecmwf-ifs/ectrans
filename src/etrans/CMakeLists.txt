

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
  set(sed_rules ${PROJECT_SOURCE_DIR}/src/etrans/sedrenames.txt)

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
  file(MAKE_DIRECTORY ${destination}/etrans_${backend})

  ecbuild_list_add_pattern( LIST absolute_files GLOB etrans/*.h SOURCE_DIR ${_PAR_INCLUDE_DIRECTORY} QUIET )
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





# TODO: move precision-independent files to common
#add_subdirectory( common )

if( HAVE_CPU)
  add_subdirectory( cpu )
endif()

# placeholder
#if( HAVE_GPU )
#  add_subdirectory( gpu )
#endif()


if (FALSE)
# original cmake file for etrans; keeping it for reference, but should be cleaned later
message(FATAL_ERROR "Hold it right there!")

# build list of sources to add to trans library 
# (using CMAKE_CURRENT_SOURCE_DIR is necessary because sources are in a different directory than the target library (trans_${prec})
ecbuild_list_add_pattern( LIST etrans_src
					  GLOB
							${CMAKE_CURRENT_SOURCE_DIR}/biper/internal/*
							${CMAKE_CURRENT_SOURCE_DIR}/biper/external/*
							${CMAKE_CURRENT_SOURCE_DIR}/etrans/internal/*
							${CMAKE_CURRENT_SOURCE_DIR}/etrans/external/*
					  QUIET
					)

# dummies to be able to loop over precisions
set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

# loop over precisions
foreach( prec sp dp )
  if( HAVE_${prec} )
	# add sources
	target_sources(trans_${prec} PRIVATE ${etrans_src})
	# add include directories
	target_include_directories(trans_${prec}
	    PUBLIC
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/biper/include>
			$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/etrans/include>
    )
  endif()
endforeach()

# install headers
file( GLOB etrans_interface biper/include/* etrans/include/*)
install(
  FILES        ${etrans_interface}
  DESTINATION  include/ectrans
)

endif()