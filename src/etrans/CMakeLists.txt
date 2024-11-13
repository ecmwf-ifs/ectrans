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
