# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Preprocess module file containing version information
configure_file( internal/common/ectrans_version_mod.F90.in internal/common/ectrans_version_mod.F90 )

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

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                algor/*
                                internal/source/*
				external/spec_norm.F90
				external/dist_grid.F90
                          QUIET
                        )
ecbuild_list_add_pattern( LIST trans_src_common
                          GLOB
                                sharedmem/*
                                internal/common/*
                          QUIET
                        )
ecbuild_list_add_pattern( LIST trans_src_ext
                          GLOB
				external/dir__trans.F90
				external/inv__trans.F90
				external/dir_transad.F90
				external/inv_transad.F90
				external/trans__inq.F90
				external/setup__trans.F90
				external/setup__trans0.F90
				external/dist_spec.F90
				external/gath_spec.F90
                                module/*
                          QUIET
                        )


function(DEFINE_CODE_SPEC) 
set (oneValueArgs INDIR OUTDIR SUBDIR BUILDTYPE SUFFIX)
set (multiValueArgs SRCFILES)
cmake_parse_arguments(DPSPREFIX "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
set(indir   ${DPSPREFIX_INDIR})
set(outdir  ${DPSPREFIX_OUTDIR})
set(subdir  ${DPSPREFIX_SUBDIR})
set(destdir ${DPSPREFIX_BUILDTYPE})
set(files   ${DPSPREFIX_SRCFILES})
set(suffix  ${DPSPREFIX_SUFFIX})

set(variantsed "s/VARIANTDESIGNATOR/"${DPSPREFIX_BUILDTYPE}"/g")
list(JOIN variantsed "" varianttmp)
separate_arguments(vartypesed  UNIX_COMMAND "${varianttmp}")

set(sourcesdir "${indir}/${subdir}")
set(destdir  "${outdir}/${destdir}/${subdir}")
message(STATUS "Sources dir:" ${sourcesdir})
message(STATUS "Destination dir:" ${destdir})
foreach( file_i ${files})
        message (STATUS "Define preprocessing for file: ${destdir}/${file_i}")
        add_custom_command( 
        COMMAND  cat ${PROJECT_SOURCE_DIR}/src/trans/include/ectrans/sedrenames.txt | 
	         sed -e ${vartypesed} | 
		 sed -f -  ${sourcesdir}/${file_i} > "${destdir}/${file_i}${suffix}"
        DEPENDS  ${sourcesdir}/${file_i}
        OUTPUT  ${destdir}/${file_i}${suffix}
        COMMENT "Preprocessing file: ${sourcesdir}/${file_i} to: ${destdir}/${file_i}${suffix}"
        VERBATIM
)
endforeach( file_i )
endfunction(DEFINE_CODE_SPEC)

function(GENERATE_TRANS_BUILD) 
set (oneValueArgs BUILDTYPE )
cmake_parse_arguments(GBPREFIX "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
set(buildtype ${GBPREFIX_BUILDTYPE}) #buildtype_arg
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/specialized_src/${buildtype}/trans/algor) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/specialized_src/${buildtype}/trans/internal/source) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/specialized_src/${buildtype}/trans/external) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/specialized_src/${buildtype}/trans/module) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/specialized_src/${buildtype}/trans/include/ectrans) 
define_code_spec(INDIR ${PROJECT_SOURCE_DIR}/src OUTDIR  "${CMAKE_BINARY_DIR}/specialized_src"  BUILDTYPE  ${buildtype}  SUBDIR "trans/" SRCFILES ${trans_src} SUFFIX ".pp.F90" )
set(trans_sources ${trans_src})
list(TRANSFORM trans_sources PREPEND ${CMAKE_BINARY_DIR}/specialized_src/${buildtype}/trans/)
list(TRANSFORM trans_sources APPEND .pp.F90)
define_code_spec(INDIR ${PROJECT_SOURCE_DIR}/src OUTDIR  "${CMAKE_BINARY_DIR}/specialized_src"  BUILDTYPE  ${buildtype}  SUBDIR "trans/" SRCFILES ${trans_src_ext} )
set(trans_sources_listname "trans_sources_${buildtype}")
set("${trans_sources_listname}" ${trans_sources} PARENT_SCOPE)
message(STATUS "Trans source list name:"  ${trans_sources_listname})

set(trans_sources_ext ${trans_src_ext})
list(TRANSFORM trans_sources_ext PREPEND ${CMAKE_BINARY_DIR}/specialized_src/${buildtype}/trans/)
set(trans_sources_ext_listname "trans_sources_ext_${buildtype}")
message(STATUS "Trans source list name:"  ${trans_sources_ext_listname})
set("${trans_sources_ext_listname}" ${trans_sources_ext} PARENT_SCOPE)

endfunction(GENERATE_TRANS_BUILD)



foreach( prec dp sp )
  if( HAVE_${prec} )
    if( ${prec} MATCHES "sp" )
       set(precdef JPRM)
    elseif (${prec} MATCHES "dp")
       set(precdef JPRD)
    endif()
    generate_trans_build( BUILDTYPE ${prec}  )
    if( NOT HAVE_FFTW )
      ecbuild_list_exclude_pattern( LIST trans_sources_${prec} REGEX tpm_fftw.F90 )
    endif()
    #message(STATUS "Preprocessed sources ext:" ${trans_sources_ext_${prec}})
    #message(STATUS "Preprocessed sources:" ${trans_sources_${prec}})
    ecbuild_add_library(
      TARGET           trans_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_sources_${prec}}
                       ${trans_src_common}
                       ${trans_sources_ext_${prec}}
                       ${CMAKE_CURRENT_BINARY_DIR}/internal/common/ectrans_version_mod.F90
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
                       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PRIVATE_INCLUDES ${CMAKE_BINARY_DIR}/specialized_src/${prec}/trans/external
      PRIVATE_DEFINITIONS SYMBOLSUFFIX=${prec}
                          JPRB=${precdef}
                          jprb=${precdef}
                          PARKIND1=EC_PARKIND
                          parkind1=EC_PARKIND
      PUBLIC_LIBS      fiat
    )

    ectrans_target_fortran_module_directory(
      TARGET            trans_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/trans_${prec}
      INSTALL_DIRECTORY module/trans_${prec}
    )
    target_link_libraries( trans_${prec} PUBLIC fiat parkind_${prec} )
    if( HAVE_FFTW )
      set( FFTW_LINK PRIVATE )
      if( LAPACK_${prec} MATCHES "mkl" AND NOT FFTW_LIBRARIES MATCHES "mkl" )
          ecbuild_warn( "Danger: Both MKL and FFTW are linked in trans_${prec}. "
                        "No guarantees on link order can be made for the final executable.")
          set( FFTW_LINK PUBLIC ) # Attempt anyway to give FFTW precedence
      endif()
      ecbuild_debug("target_link_libraries( trans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )")
      target_link_libraries( trans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )
      target_include_directories( trans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
      target_compile_definitions( trans_${prec} PRIVATE WITH_FFTW )
    endif()
    ecbuild_debug("target_link_libraries( trans_${prec} ${LAPACK_LINK} ${LAPACK_${prec}} )")
    target_link_libraries( trans_${prec} ${LAPACK_LINK} ${LAPACK_${prec}} )
    if( HAVE_OMP )
      ecbuild_debug("target_link_libraries( trans_${prec} PRIVATE OpenMP::OpenMP_Fortran )")
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
