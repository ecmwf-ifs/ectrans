# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Apply workarounds for some known compilers
## see trans/ for example

function(generate_backend_sources)
  set (options)
  set (oneValueArgs BACKEND DESTINATION OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})
  file(MAKE_DIRECTORY ${destination}/biper/internal)
  file(MAKE_DIRECTORY ${destination}/biper/external)
  file(MAKE_DIRECTORY ${destination}/internal)
  file(MAKE_DIRECTORY ${destination}/external)

  ecbuild_list_add_pattern( LIST files
    GLOB
    internal/*.F90
    external/*.F90
    biper/internal/*.F90
          biper/external/*.F90
    QUIET
  )

  set(outfiles)
  foreach(file_i ${files})
    get_filename_component(outfile_name    ${file_i} NAME)
    get_filename_component(outfile_name_we ${file_i} NAME_WE)
    get_filename_component(outfile_ext     ${file_i} EXT)
    get_filename_component(outfile_dir     ${file_i} DIRECTORY)
    set(outfile "${destination}/${file_i}")
    ecbuild_debug("Generate ${outfile}")
    generate_file(BACKEND ${backend} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile})
    list(APPEND outfiles ${outfile})
  endforeach(file_i)
  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)

set( BUILD_INTERFACE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include/ectrans )

foreach( prec dp sp )
  if( HAVE_${prec} )

    generate_backend_includes(BACKEND ${prec} TARGET ectrans_etrans_${prec}_includes DESTINATION ${BUILD_INTERFACE_INCLUDE_DIR} INCLUDE_DIRECTORY ${PROJECT_SOURCE_DIR}/src/etrans/include )
    generate_backend_sources( BACKEND ${prec} OUTPUT ectrans_etrans_${prec}_src  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_etrans_${prec})
    ecbuild_add_library(
      TARGET           ectrans_etrans_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${ectrans_etrans_${prec}_src}
      PUBLIC_INCLUDES  $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
                       $<BUILD_INTERFACE:${BUILD_INTERFACE_INCLUDE_DIR}>
      PUBLIC_LIBS      fiat ectrans_common ectrans_${prec}_includes ectrans_${prec} ectrans_etrans_${prec}_includes
    )

    ectrans_target_fortran_module_directory(
      TARGET            ectrans_etrans_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans
      INSTALL_DIRECTORY module/ectrans
    )

    set( FFTW_LINK PRIVATE )
    if( LAPACK_LIBRARIES MATCHES "mkl" AND NOT FFTW_LIBRARIES MATCHES "mkl" )
        ecbuild_warn( "Danger: Both MKL and FFTW are linked in trans_${prec}. "
                      "No guarantees on link order can be made for the final executable.")
        set( FFTW_LINK PUBLIC ) # Attempt anyway to give FFTW precedence
    endif()
    ecbuild_debug("target_link_libraries( trans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )")
    target_link_libraries( ectrans_etrans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )
    target_include_directories( ectrans_etrans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
    target_compile_definitions( ectrans_etrans_${prec} PRIVATE WITH_FFTW )
    # daand: lam transforms don't need lapack
    #ecbuild_debug("target_link_libraries( ectrans_etrans_${prec} PRIVATE ${LAPACK_LIBRARIES} )")
    #target_link_libraries( ectrans_${prec} PRIVATE ${LAPACK_LIBRARIES} )

    if( HAVE_OMP )
      ecbuild_debug("target_link_libraries( ectrans_${prec} PRIVATE OpenMP::OpenMP_Fortran )")
      target_link_libraries( ectrans_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    # This interface library is for backward compatibility, and provides the older includes
    ecbuild_add_library( TARGET etrans_${prec} TYPE INTERFACE )
    target_include_directories( etrans_${prec} INTERFACE $<BUILD_INTERFACE:${BUILD_INTERFACE_INCLUDE_DIR}/etrans_${prec}> )
    target_include_directories( etrans_${prec} INTERFACE $<INSTALL_INTERFACE:include/ectrans/etrans_${prec}> )
    target_link_libraries( trans_${prec} INTERFACE fiat ectrans_${prec} ectrans_etrans_${prec} parkind_${prec})
  endif()
endforeach()

## Install trans interface
install( DIRECTORY ${BUILD_INTERFACE_INCLUDE_DIR}/ DESTINATION include/ectrans )
