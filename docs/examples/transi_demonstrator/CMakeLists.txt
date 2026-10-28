cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( ectrans_demonstrator LANGUAGES C )

find_package( fiat REQUIRED )
find_package( ectrans REQUIRED )

add_executable( transi_test transi_test.c )
target_link_libraries( transi_test PRIVATE transi_dp )

