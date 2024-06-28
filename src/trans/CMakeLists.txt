# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Preprocess module file containing version information
configure_file( internal/ectrans_version_mod.F90.in internal/ectrans_version_mod.F90 )

## Apply workarounds for some known compilers

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8.7 )

    # Fix for IFS "CONGRAD: SPTSV/DPTSV returned non-zero info with crayftn 8.7.7 (cdt/18.12)
    ectrans_add_compile_options(
        SOURCES internal/ftinv_ctlad_mod.F90
        FLAGS   "-O0,fp1,omp")

  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 18 AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19 )

    # See https://github.com/ecmwf-ifs/ectrans/issues/17
    ectrans_add_compile_options(
        SOURCES algor/butterfly_alg_mod.F90
        FLAGS   "-check nopointers")
  endif()
endif()

## Sources which are precision independent can go into a common library
list( APPEND ectrans_common_src
  algor/ectrans_blas_mod.F90
  sharedmem/sharedmem_mod.F90
  sharedmem/sharedmem.c
  internal/abort_trans_mod.F90
  internal/cpledn_mod.F90
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
  ${CMAKE_CURRENT_BINARY_DIR}/internal/ectrans_version_mod.F90
)
list( APPEND ectrans_common_includes
  include/ectrans/get_current.h
  include/ectrans/setup_trans0.h
)

ecbuild_add_library(
  TARGET           ectrans_common
  LINKER_LANGUAGE  Fortran
  SOURCES          ${ectrans_common_src}
  PUBLIC_LIBS      fiat
  PRIVATE_LIBS     ${LAPACK_LIBRARIES}
)
ectrans_target_fortran_module_directory(
  TARGET            ectrans_common
  MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans
  INSTALL_DIRECTORY module/ectrans
)

if( HAVE_OMP )
  ecbuild_debug("target_link_libraries( trans_${prec} PRIVATE OpenMP::OpenMP_Fortran )")
  target_link_libraries( ectrans_common PRIVATE OpenMP::OpenMP_Fortran )
endif()



function(generate_file)
  set (options)
  set (oneValueArgs INPUT OUTPUT BACKEND)
  set (multiValueArgs)
  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(output  ${_PAR_OUTPUT})
  set(input   ${_PAR_INPUT})
  set(backend ${_PAR_BACKEND})

  set( JPRB_dp JPRD )
  set( jprb_dp jprd )
  set( JPRB_sp JPRM )
  set( jprb_sp jprm )

  add_custom_command(
    OUTPUT  ${output}
    COMMAND  cat ${CMAKE_CURRENT_SOURCE_DIR}/sedrenames.txt |
             sed -e "s/VARIANTDESIGNATOR/${backend}/g" |
             sed -e "s/TYPEDESIGNATOR_UPPER/${JPRB_${backend}}/g" |
             sed -e "s/TYPEDESIGNATOR_LOWER/${jprb_${backend}}/g" |
             sed -rf - ${CMAKE_CURRENT_SOURCE_DIR}/${input} > ${output}
    DEPENDS  ${input} ${CMAKE_CURRENT_SOURCE_DIR}/sedrenames.txt
    COMMENT "Generating ${output}"
    VERBATIM
  )
endfunction(generate_file)


function(generate_backend_includes)
  set (options)
  set (oneValueArgs BACKEND TARGET DESTINATION)
  set (multiValueArgs)
  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(destination ${_PAR_DESTINATION} )
  set(backend     ${_PAR_BACKEND})

  file(MAKE_DIRECTORY ${destination})
  file(MAKE_DIRECTORY ${destination}/trans_${backend})

  ecbuild_list_add_pattern( LIST files GLOB include/ectrans/*.h QUIET )

  set( outfiles )
  foreach(file_i ${files})
    get_filename_component(outfile_name    ${file_i} NAME)
    get_filename_component(outfile_name_we ${file_i} NAME_WE)
    get_filename_component(outfile_ext     ${file_i} EXT)
    get_filename_component(outfile_dir     ${file_i} DIRECTORY)
    if (${file_i} IN_LIST ectrans_common_includes)
      configure_file(${file_i} ${destination}/${outfile_name})
    else()
      set(outfile "${destination}/${outfile_name_we}_${backend}${outfile_ext}")
      ecbuild_debug("Generate ${outfile}")
      generate_file(BACKEND ${backend} INPUT ${file_i} OUTPUT ${outfile})
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



function(generate_backend_sources)
  set (options)
  set (oneValueArgs BACKEND DESTINATION OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})
  file(MAKE_DIRECTORY ${destination}/algor)
  file(MAKE_DIRECTORY ${destination}/internal)
  file(MAKE_DIRECTORY ${destination}/external)

  ecbuild_list_add_pattern( LIST files
    GLOB
          algor/*
          internal/*
          external/*
    QUIET
  )

  set(outfiles)
  foreach(file_i ${files})
    if(NOT (${file_i} IN_LIST ectrans_common_src))
      get_filename_component(outfile_name    ${file_i} NAME)
      get_filename_component(outfile_name_we ${file_i} NAME_WE)
      get_filename_component(outfile_ext     ${file_i} EXT)
      get_filename_component(outfile_dir     ${file_i} DIRECTORY)
      set(outfile "${destination}/${file_i}")
      ecbuild_debug("Generate ${outfile}")
      generate_file(BACKEND ${backend} INPUT ${file_i} OUTPUT ${outfile})
      list(APPEND outfiles ${outfile})
    endif()
  endforeach(file_i)
  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)

set( BUILD_INTERFACE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include/ectrans )

foreach( prec dp sp )
  if( HAVE_${prec} )

    generate_backend_includes(BACKEND ${prec} TARGET ectrans_${prec}_includes DESTINATION ${BUILD_INTERFACE_INCLUDE_DIR} )
    generate_backend_sources( BACKEND ${prec} OUTPUT ectrans_${prec}_src  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_${prec})
    ecbuild_add_library(
      TARGET           ectrans_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${ectrans_${prec}_src}
      PUBLIC_INCLUDES  $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
                       $<BUILD_INTERFACE:${BUILD_INTERFACE_INCLUDE_DIR}>
      PUBLIC_LIBS      ectrans_common ectrans_${prec}_includes
    )

    ectrans_target_fortran_module_directory(
      TARGET            ectrans_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans_${prec}
      INSTALL_DIRECTORY module/ectrans_${prec}
    )
    target_link_libraries( ectrans_${prec} PUBLIC fiat)

    set( FFTW_LINK PRIVATE )
    if( LAPACK_LIBRARIES MATCHES "mkl" AND NOT FFTW_LIBRARIES MATCHES "mkl" )
        ecbuild_warn( "Danger: Both MKL and FFTW are linked in trans_${prec}. "
                      "No guarantees on link order can be made for the final executable.")
        set( FFTW_LINK PUBLIC ) # Attempt anyway to give FFTW precedence
    endif()
    ecbuild_debug("target_link_libraries( trans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )")
    target_link_libraries( ectrans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )
    target_include_directories( ectrans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
    target_compile_definitions( ectrans_${prec} PRIVATE WITH_FFTW )
    ecbuild_debug("target_link_libraries( ectrans_${prec} PRIVATE ${LAPACK_LIBRARIES} )")
    target_link_libraries( ectrans_${prec} PRIVATE ${LAPACK_LIBRARIES} )

    if( HAVE_OMP )
      ecbuild_debug("target_link_libraries( ectrans_${prec} PRIVATE OpenMP::OpenMP_Fortran )")
      target_link_libraries( ectrans_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    # This interface library is for backward compatibility, and provides the older includes
    ecbuild_add_library( TARGET trans_${prec} TYPE INTERFACE )
    target_include_directories( trans_${prec} INTERFACE $<BUILD_INTERFACE:${BUILD_INTERFACE_INCLUDE_DIR}/trans_${prec}> )
    target_include_directories( trans_${prec} INTERFACE $<INSTALL_INTERFACE:include/ectrans/trans_${prec}> )
    target_link_libraries( trans_${prec} INTERFACE fiat ectrans_${prec} parkind_${prec})
  endif()
endforeach()

## Install trans interface
install( DIRECTORY ${BUILD_INTERFACE_INCLUDE_DIR}/ DESTINATION include/ectrans )
