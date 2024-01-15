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

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                sharedmem/*
                                algor/*
                                internal/*
				#external/spec_norm.F90
                                module/*
				#  ${CMAKE_CURRENT_BINARY_DIR}/internal/ectrans_version_mod.F90
                          QUIET
                        )
ecbuild_list_add_pattern( LIST trans_src_ext
                          GLOB
				external/spec_norm.F90
				external/dir__trans.F90
				external/inv__trans.F90
				external/trans__inq.F90
				external/setup__trans.F90
				external/setup__trans0.F90
                          QUIET
                        )

get_target_property(fiat_include fiat INTERFACE_INCLUDE_DIRECTORIES)
set(fiat_include_opts ${fiat_include})
list(TRANSFORM fiat_include_opts PREPEND " -I")
list(JOIN fiat_include_opts " " fiat_include_opts_j1)
separate_arguments(fiat_include_opts_j  UNIX_COMMAND "${fiat_include_opts_j1}")


function(DEFINE_PREPROCESS) 
set (options ONGPU )
set (oneValueArgs INDIR OUTDIR SUBDIR BUILDTYPE PRECTYPE)
set (multiValueArgs SOURCEFILES SPECOPTS)
cmake_parse_arguments(DPSPREFIX "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
set(indir   ${DPSPREFIX_INDIR})
set(outdir  ${DPSPREFIX_OUTDIR})
set(subdir  ${DPSPREFIX_SUBDIR})
set(destdir ${DPSPREFIX_BUILDTYPE})
set(custopt ${DPSPREFIX_SPECOPTS})
set(files   ${DPSPREFIX_SOURCEFILES})

set(precsed "s/JPRBC/"${DPSPREFIX_PRECTYPE}"/g")
list(JOIN precsed "" precsedtmp)
separate_arguments(prectypesed  UNIX_COMMAND "${precsedtmp}")

set(variantsed "s/VARIANTDESIGNATOR/"${DPSPREFIX_BUILDTYPE}"/g")
list(JOIN variantsed "" varianttmp)
separate_arguments(vartypesed  UNIX_COMMAND "${varianttmp}")

message (STATUS "Subdir: " ${subdir})
message (STATUS "Buildtype: " ${destdir})
message (STATUS "Specopts: " ${custopt})
message (STATUS "Files " "${files}")
set(sourcesdir "${indir}/${subdir}")
set(destdir  "${outdir}/${destdir}/${subdir}")
message(STATUS "Sources dir:" ${sourcesdir})
message(STATUS "Destination dir:" ${destdir})
foreach( file_i ${files})
        message (STATUS "Define preprocessing for file: ${destdir}/${file_i}")
        add_custom_command( 
#	COMMAND gcc-13 -E -P  -cpp  ${sourcesdir}/${file_i} -I${FFTW_INCLUDE_DIRS} ${fiat_include_opts_j}  -I${PROJECT_SOURCE_DIR}/src/trans/include/ectrans -I${PROJECT_SOURCE_DIR}/src/trans/external -DSYMBOLSUFFIX=${prec} -DPRECOPT=${precno} > "${destdir}/${file_i}.pp.f90"
COMMAND  cat ${PROJECT_SOURCE_DIR}/src/trans/include/ectrans/sedrenames.txt | sed -e ${vartypesed} -e ${prectypesed}  > sedtest.txt && sed -f sedtest.txt  ${sourcesdir}/${file_i} > "${destdir}/${file_i}.pp.F90"

        DEPENDS  ${sourcesdir}/${file_i}
        OUTPUT  ${destdir}/${file_i}.pp.F90
        COMMENT "Preprocessing ${sourcesdir}/${file_i} to ${destdir}/${file_i}.pp.F90"
        VERBATIM
)
endforeach( file_i )
endfunction(DEFINE_PREPROCESS)
function(DEFINE_PREPROCESS_EXT) 
set (options ONGPU )
set (oneValueArgs INDIR OUTDIR SUBDIR BUILDTYPE PRECTYPE)
set (multiValueArgs SOURCEFILES SPECOPTS)
cmake_parse_arguments(DPSPREFIX "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
set(indir   ${DPSPREFIX_INDIR})
set(outdir  ${DPSPREFIX_OUTDIR})
set(subdir  ${DPSPREFIX_SUBDIR})
set(destdir ${DPSPREFIX_BUILDTYPE})
set(custopt ${DPSPREFIX_SPECOPTS})
set(files   ${DPSPREFIX_SOURCEFILES})

set(precsed "s/JPRBC/"${DPSPREFIX_PRECTYPE}"/g")
list(JOIN precsed "" precsedtmp)
separate_arguments(prectypesed  UNIX_COMMAND "${precsedtmp}")

set(variantsed "s/VARIANTDESIGNATOR/"${DPSPREFIX_BUILDTYPE}"/g")
list(JOIN variantsed "" varianttmp)
separate_arguments(vartypesed  UNIX_COMMAND "${varianttmp}")

message (STATUS "Subdir: " ${subdir})
message (STATUS "Buildtype: " ${destdir})
message (STATUS "Specopts: " ${custopt})
message (STATUS "Files " "${files}")
set(sourcesdir "${indir}/${subdir}")
set(destdir  "${outdir}/${destdir}/${subdir}")
message(STATUS "Sources dir:" ${sourcesdir})
message(STATUS "Destination dir:" ${destdir})
foreach( file_i ${files})
        message (STATUS "Define preprocessing for file: ${destdir}/${file_i}")
        add_custom_command( 
#	COMMAND gcc-13 -E -P  -cpp  ${sourcesdir}/${file_i} -I${FFTW_INCLUDE_DIRS} ${fiat_include_opts_j}  -I${PROJECT_SOURCE_DIR}/src/trans/include/ectrans -I${PROJECT_SOURCE_DIR}/src/trans/external -DSYMBOLSUFFIX=${prec} -DPRECOPT=${precno} > "${destdir}/${file_i}.pp.f90"
COMMAND  cat ${PROJECT_SOURCE_DIR}/src/trans/include/ectrans/sedrenamesext.txt | sed -e ${vartypesed} -e ${prectypesed}  > sedtest.txt && sed -f sedtest.txt  ${sourcesdir}/${file_i} > "${destdir}/${file_i}"

        DEPENDS  ${sourcesdir}/${file_i}
        OUTPUT  ${destdir}/${file_i}
        COMMENT "Preprocessing ${sourcesdir}/${file_i} to ${destdir}/${file_i}"
        VERBATIM
)
endforeach( file_i )
endfunction(DEFINE_PREPROCESS_EXT)

function(GENERATE_BUILD) 
set (options ONGPU )
set (oneValueArgs BUILDTYPE PRECTYPE)
set (multiValueArgs SPECOPTS )
cmake_parse_arguments(GBPREFIX "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
set(buildtype ${GBPREFIX_BUILDTYPE}) #buildtype_arg
set(specopts  ${GBPREFIX_SPECOPTS}) #specopts_arg} )
set(ongpu ${GBPREFIX_ONGPU}) #  ${driverfile_arg} )
set(precd ${GBPREFIX_PRECTYPE})
message (STATUS "Buildtype: " ${buildtype})
message (STATUS "Specopts: " ${specopts})
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src) 
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}) 
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/sharedmem) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/algor) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/internal) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/external) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/module) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/include/ectrans) 
list(REMOVE_ITEM trans_src internal/ectrans_version_mod.F90.in) 
list(REMOVE_ITEM trans_src internal/ectrans_version_mod.F90) 
list(REMOVE_ITEM trans_src sharedmem/sharedmem.c) 
define_preprocess(INDIR ${PROJECT_SOURCE_DIR} OUTDIR  "${CMAKE_BINARY_DIR}/preprocessed_src"  SUBDIR "src/trans/"  BUILDTYPE  ${buildtype} PRECTYPE ${precd} SOURCEFILES ${trans_src}  SPECOPTS ${specopts})
define_preprocess_ext(INDIR ${PROJECT_SOURCE_DIR} OUTDIR  "${CMAKE_BINARY_DIR}/preprocessed_src"  SUBDIR "src/trans/"  BUILDTYPE  ${buildtype} PRECTYPE ${precd} SOURCEFILES ${trans_src_ext}  SPECOPTS ${specopts})
set(trans_sources ${trans_src})
set(trans_sources_ext ${trans_src_ext})
#list(REMOVE_ITEM trans_sources ${CMAKE_CURRENT_BINARY_DIR}/internal/ectrans_version_mod.F90) 
#list(REMOVE_ITEM trans_sources sharedmem/sharedmem.c) 
list(TRANSFORM trans_sources PREPEND ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/)
list(TRANSFORM trans_sources APPEND .pp.F90)
list(TRANSFORM trans_sources_ext PREPEND ${CMAKE_BINARY_DIR}/preprocessed_src/${buildtype}/src/trans/)
list(APPEND trans_sources internal/ectrans_version_mod.F90) 
list(APPEND trans_sources sharedmem/sharedmem.c) 
set(trans_sources_listname "trans_sources_${buildtype}")
message(STATUS "Trans source list name:"  ${trans_sources_listname})
set("${trans_sources_listname}" ${trans_sources} PARENT_SCOPE)
set(trans_sources_ext_listname "trans_sources_ext_${buildtype}")
message(STATUS "Trans source list name:"  ${trans_sources_ext_listname})
set("${trans_sources_ext_listname}" ${trans_sources_ext} PARENT_SCOPE)
endfunction(GENERATE_BUILD)




foreach( prec dp sp )
    if( ${prec} MATCHES "sp" )
       set(precno 1)
    elseif (${prec} MATCHES "dp")
       set(precno 2)
    endif()
endforeach( prec )


if( NOT HAVE_FFTW )
  ecbuild_list_exclude_pattern( LIST trans_src REGEX tpm_fftw.F90 )
endif()

foreach( prec dp sp )
  if( HAVE_${prec} )
  
    if( ${prec} MATCHES "sp" )
       set(precno 1)
       set(precdef 'JPRM')
    elseif (${prec} MATCHES "dp")
       set(precno 2)
       set(precdef 'JPRD')
    endif()
    generate_build( SPECOPTS  ${custopts} BUILDTYPE ${prec} PRECTYPE ${precdef} )
    message(STATUS "Preprocessed sources:" ${trans_sources_${prec}})
    message(STATUS "Preprocessed sources ext:" ${trans_sources_ext_${prec}})

    ecbuild_add_library(
      TARGET           trans_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_sources_${prec}}
                       ${trans_sources_ext_${prec}}
                       ${CMAKE_CURRENT_BINARY_DIR}/internal/ectrans_version_mod.F90
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include>
                       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PRIVATE_INCLUDES ${CMAKE_BINARY_DIR}/preprocessed_src/${prec}/src/trans/external
      PRIVATE_DEFINITIONS SYMBOLSUFFIX=${prec}
                               PRECOPT=${precno}
      PUBLIC_LIBS      fiat parkind_${prec} 
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
