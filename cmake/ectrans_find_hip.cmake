# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ectrans_find_hip )
  # This macro finds all HIP related libraries, if found, HAVE_HIP=1

  set( options "" )
  set( single_value_args REQUIRED )
  set( multi_value_args "" )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  set(HIP_REQUIRED "")
  if( _PAR_REQUIRED )
    set(HIP_REQUIRED "REQUIRED" )
  endif()

  set(HAVE_HIP 1)

  # Setup ROCM_PATH
  if (NOT DEFINED ROCM_PATH )
    find_path(ROCM_PATH
      hip
      ENV{ROCM_DIR}
      ENV{ROCM_PATH}
      ENV{HIP_PATH}
      ${HIP_PATH}/..
      ${HIP_ROOT_DIR}/../
      ${ROCM_ROOT_DIR}
      /opt/rocm)
  endif()
  ecbuild_info("ROCM path: ${ROCM_PATH}")
  # Update CMAKE_PREFIX_PATH to make sure all the configs that hip depends on are found.
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${ROCM_PATH}")

  set(HAVE_HIP 1)

  ecbuild_add_option( FEATURE HIP_LANGUAGE DEFAULT ON )
  
  if( HAVE_HIP_LANGUAGE )
    include(CheckLanguage)
    check_language(HIP)
    if(NOT CMAKE_HIP_COMPILER)
      if( _PAR_REQUIRED )
        ecbuild_error("HIP compiler not found")
      else()
        ecbuild_info("HIP compiler not found: HAVE_HIP=0")
        set(HAVE_HIP 0)
      endif()
    else()
      enable_language(HIP)
      ecbuild_info("HIP compiler found: ${CMAKE_HIP_COMPILER}")
      ecbuild_info("HIP target architecture: ${CMAKE_HIP_ARCHITECTURES}")
    endif()

    # Find HIP libraries
    find_package(hip REQUIRED CONFIG)
    if( NOT hip_FOUND )
      ecbuild_info("hip libraries not found: HAVE_HIP=0")
      set( HAVE_HIP 0 )
    endif()

    ecbuild_info("HIP version: ${hip_VERSION}")

  else()
    ecbuild_info("HIP sources will be compiled with C++ compiler with added flags")
    ecbuild_info("HIP target architecture: ${CMAKE_HIP_ARCHITECTURES}")
    enable_language(CXX)
    set(CMAKE_MODULE_PATH $ENV{HIP_ROOT}/cmake ${CMAKE_MODULE_PATH})
    find_package(HIP)
    if ( NOT HIP_FOUND )
      ecbuild_info("HIP not found: HAVE_HIP=0")
      set( HAVE_HIP 0)
    endif()
    ecbuild_info("HIP version: ${HIP_VERSION}")

  endif()

  if( HAVE_HIP )
    # Find HIP libraries
    # find_package(hip REQUIRED CONFIG)
    # if( NOT hip_FOUND )
    #   ecbuild_info("hip libraries not found: HAVE_HIP=0")
    #   set( HAVE_HIP 0 )
    # endif()

    find_package(hipblas CONFIG ${HIP_REQUIRED})
    if( NOT hipblas_FOUND )
      ecbuild_info("hipblas libraries not found: HAVE_HIP=0")
      set( HAVE_HIP 0 )
    endif()

    find_package(hipfft  CONFIG ${HIP_REQUIRED})
    if( NOT hipfft_FOUND )
      ecbuild_info("hipfft libraries not found: HAVE_HIP=0")
      set( HAVE_HIP 0 )
    endif()

    find_package(rocblas CONFIG ${HIP_REQUIRED})
    if( NOT rocblas_FOUND )
      ecbuild_info("rocblas libraries not found: HAVE_HIP=0")
      set( HAVE_HIP 0 )
    endif()

    find_package(rocfft  CONFIG ${HIP_REQUIRED})
    if( NOT rocfft_FOUND )
      ecbuild_info("rocfft libraries not found: HAVE_HIP=0")
      set( HAVE_HIP 0 )
    endif()

    if( HAVE_HIP )
      list( APPEND ECTRANS_GPU_HIP_LIBRARIES ${hipblas_LIBRARIES} ${hipfft_LIBRARIES})
      list( APPEND ECTRANS_GPU_HIP_LIBRARIES ${rocblas_LIBRARIES} ${rocfft_LIBRARIES})
    endif()

  endif()
  ecbuild_info("HIP libraries: ${ECTRANS_GPU_HIP_LIBRARIES}")
endmacro()

macro( ectrans_declare_hip_sources )
  set( options QUIET )
  set( single_value_args "" )
  set( multi_value_args SOURCES SOURCES_GLOB )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )
  set( source_files ${_PAR_SOURCES} )
  if( _PAR_SOURCES_GLOB )
    ecbuild_list_add_pattern( LIST source_files
                              GLOB ${_PAR_SOURCES_GLOB}
                              QUIET )
  endif()

  if( HAVE_HIP_LANGUAGE )
    if(NOT _PAR_QUIET)
      ecbuild_info("Applying HIP language to ${source_files}")
    endif()
    set_source_files_properties( ${source_files} PROPERTIES LANGUAGE HIP )
  else()
    if(NOT _PAR_QUIET)
      ecbuild_info("Applying HIP flags to ${source_files}")
    endif()
    set_source_files_properties( ${source_files} PROPERTIES LANGUAGE CXX )
    set( _flags "-x hip" )
    if( CMAKE_HIP_FLAGS )
      set( _flags "${_flags} ${CMAKE_HIP_FLAGS}" )
    endif()
    if( CMAKE_HIP_ARCHITECTURES )
      set( _flags "${_flags} --offload-arch=${CMAKE_HIP_ARCHITECTURES}" )
    endif()
    set_source_files_properties( ${source_files} PROPERTIES COMPILE_FLAGS "${_flags}" )
  endif()
endmacro()
