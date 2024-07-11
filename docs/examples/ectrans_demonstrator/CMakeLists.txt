# (C) Copyright 2024- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( ectrans_demonstrator LANGUAGES Fortran )

find_package( fiat REQUIRED )
find_package( ectrans REQUIRED )

add_executable( ectrans_demonstrator ectrans_demonstrator.F90 )

target_include_directories( ectrans_demonstrator 
                            PRIVATE ${CMAKE_SOURCE_DIR}/ectrans/src/trans/include/ectrans )
target_include_directories( ectrans_demonstrator PRIVATE ${fiat_ROOT}/module/parkind_sp )
target_link_libraries( ectrans_demonstrator PRIVATE trans_sp)
