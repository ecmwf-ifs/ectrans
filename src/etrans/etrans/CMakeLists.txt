
## Assemble sources

ecbuild_list_add_pattern( LIST etrans_src
                          GLOB
                                internal/*
                                external/*
								aux/*.F90
                          QUIET
                        )

set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

foreach( prec sp dp )
  if( HAVE_${prec} )

    ecbuild_add_library(
      TARGET           etrans_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${etrans_src}
      PUBLIC_INCLUDES  #$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/etrans>
                       $<INSTALL_INTERFACE:include/ectrans>
                       #$<INSTALL_INTERFACE:include>
      PUBLIC_LIBS      fiat parkind_${prec}
	  PRIVATE_LIBS     trans_${prec} biper_${prec}
    )
    ectrans_target_fortran_module_directory(
      TARGET            etrans_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/etrans_${prec}
      INSTALL_DIRECTORY module/etrans_${prec}
    )
	#target_link_libraries( biper_${prec} PUBLIC fiat parkind_${prec} trans_${prec} biper_${prec})
    #if( HAVE_FFTW )   # already resolved from trans, I presume
    #  target_link_libraries( etrans_${prec} ${FFTW_LINK} ${FFTW_LIBRARIES} )
    #  target_include_directories( etrans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
    #  target_compile_definitions( etrans_${prec} PRIVATE WITH_FFTW )
    #endif()
    #target_link_libraries( etrans_${prec} ${LAPACK_LINK} ${LAPACK_${prec}} )    # lapack isn't used by etrans
    #if( HAVE_OMP )   # already resolved from trans, I presume
    #  target_link_libraries( etrans_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    #endif()

  endif()
endforeach()

## Install trans interface

file( GLOB etrans_interface include/etrans/* )
install(
  FILES        ${etrans_interface}
  DESTINATION  include/ectrans
)
