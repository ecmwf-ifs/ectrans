macro( ectrans_find_cuda )
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      ecbuild_info("CMAKE_CUDA_ARCHITECTURES not defined, using 80")
      set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
    check_language(CUDA)
    if ( NOT CMAKE_CUDA_COMPILER )
      set( HAVE_CUDA 0 )
    else()
      enable_language(CUDA)
      set( HAVE_CUDA 1 )
      find_package( CUDAToolkit )
      if( NOT TARGET CUDA::cublas )
        ecbuild_info("No target CUDA::cublas")
        set( HAVE_CUDA 0 )
      endif()
      if( NOT TARGET CUDA::cufft )
        ecbuild_info("No target CUDA::cufft")
        set( HAVE_CUDA 1 )
      endif()
      ecbuild_info( "cuda arch               : [${CMAKE_CUDA_ARCHITECTURES}]" )
      ecbuild_info( "cublas                  : [${CUDA_cublas_LIBRARY}]" )
      ecbuild_info( "cufft                   : [${CUDA_cufft_LIBRARY}]" )

      if( TARGET CUDA::cublas )
         # When looking for CUDAToolkit, the various components are found using find_library
         # which sets the target to the full path, but the parent directory is not added
         # to the target's LINK_DIRECTORIES property. Whilst in principle this should be enough
         # to link to the CUDA libraries, as the full path appears on the link line, it seems
         # there is a transient dependency somewhere that also requires the directory path to be set.
         # This symptom only manifests on systems where libcudart lives in a separate location to
         # libcublas. It is unclear whether this is an nvhpc or a cmake bug, but nevertheless manually
         # specifying the link directory here should be harmless.
         cmake_path(GET CUDA_cublas_LIBRARY PARENT_PATH _cublas_path)
         target_link_directories( CUDA::cublas INTERFACE ${_cublas_path} )
      endif()
    endif()
endmacro()
