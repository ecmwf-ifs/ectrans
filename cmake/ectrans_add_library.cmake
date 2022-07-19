function( ectrans_add_library )

    set( options )
    set( single_value_args TYPE TARGET )
    set( multi_value_args SOURCES LINKER_LANGUAGE )
    cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN} )


    if( NOT DEFINED _PAR_SOURCES )
        # ecbuild_add_library requires a "SOURCES" keyword, which should not be required for linking an object library
        if( NOT DEFINED _PAR_LINKER_LANGUAGE )
            ecbuild_error("LINKER_LANGUAGE keyword missing")
        endif()
        set( dummy_source ${CMAKE_CURRENT_BINARY_DIR}/${_PAR_TARGET}_dummy.F90 )
        if( NOT EXISTS ${dummy_source} )
            file( WRITE ${dummy_source} "subroutine ${_PAR_TARGET}_dummy; end subroutine" )
        endif()
        ecbuild_add_library( ${ARGV} SOURCES ${dummy_source} )
    else()
        ecbuild_add_library( ${ARGV} )
    endif()

    if( _PAR_TYPE STREQUAL "OBJECT" )
        # ecbuild support for exporting object libraries is missing. For now do it manually.
        install( TARGETS ${_PAR_TARGET}
                 EXPORT  ${PROJECT_NAME}-targets
                 RUNTIME DESTINATION ${INSTALL_BIN_DIR}
                 LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
                 ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
        )
        export( TARGETS ${_PAR_TARGET} APPEND FILE "${PROJECT_TARGETS_FILE}" )
    endif()

endfunction()

