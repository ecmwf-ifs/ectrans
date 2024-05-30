## Assemble sources
ecbuild_list_add_pattern( LIST biper_src
                          GLOB
                                internal/*
                                external/*
                          QUIET
                        )

set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

foreach( prec sp dp )
  if( HAVE_${prec} )

    ecbuild_add_library(
      TARGET           biper_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${biper_src}
      #PUBLIC_INCLUDES  #$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                       #$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/biper>
                       #$<INSTALL_INTERFACE:include/biper>
                       #$<INSTALL_INTERFACE:include>
      PUBLIC_LIBS      fiat parkind_${prec}
	  PRIVATE_LIBS     trans_${prec}
    )
	
	#target_link_libraries( biper_${prec} PUBLIC fiat parkind_${prec} trans_${prec})
	
	# not sure if modules should be installed: shouldn't biper be accessed through interface routines?
    ectrans_target_fortran_module_directory(
      TARGET            biper_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/biper_${prec}
      INSTALL_DIRECTORY module/biper_${prec}
    )

  endif()
endforeach()

## Install biper interface
file( GLOB biper_interface include/biper/* )
install(
  FILES        ${biper_interface}
  DESTINATION  include/ectrans
)
