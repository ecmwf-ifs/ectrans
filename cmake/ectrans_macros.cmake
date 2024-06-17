# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( CMAKE_VERSION VERSION_LESS 3.22.0 )
    # FindCUDAToolkit from cmake < 3.22 does not support recent NVHPC directory structures
    set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/FindCUDAToolkit-cmake-3.24 ${CMAKE_MODULE_PATH} )
endif()

### Workaround to extract GIT_SHA1 from parent directory
if( NOT ${PROJECT_NAME}_GIT_SHA1 )
    get_filename_component( PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY )
    if( EXISTS ${PARENT_DIR}/.git )
        get_filename_component( PARENT_REPOSITORY_NAME ${PARENT_DIR} NAME_WE )
        get_git_head_revision( GIT_REFSPEC ${PROJECT_NAME}_GIT_SHA1 )
        string( SUBSTRING "${${PROJECT_NAME}_GIT_SHA1}" 0 7 ${PROJECT_NAME}_GIT_SHA1_SHORT )
        set( ${PROJECT_NAME}_GIT_SHA1_SHORT "${PARENT_REPOSITORY_NAME}/${${PROJECT_NAME}_GIT_SHA1_SHORT}" )
        set( ${PROJECT_NAME}_GIT_SHA1       "${PARENT_REPOSITORY_NAME}/${${PROJECT_NAME}_GIT_SHA1}" )
    endif()
endif()

include( ectrans_target_fortran_module_directory )
include( ectrans_find_lapack )
include( ectrans_find_cuda )
include( ectrans_find_hip )
include( CheckLanguage )

