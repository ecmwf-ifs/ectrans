# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Tests for analytic correctness of INV_TRANS and DIR_TRANS

foreach( precision ${precisions} )

  foreach( platform ${platforms} )

    if( "${platform}" MATCHES "cpu" )
      set( platform_tag "" )
    endif()
    if( "${platform}" MATCHES "gpu" )
      set( platform_tag "_gpu" )
    endif()

    ecbuild_add_executable(
      TARGET ectrans_test_analytic_${platform}_${precision}
      SOURCES analytic_test.F90
      LIBS fiat trans${platform_tag}_${precision} parkind_${precision}
      LINKER_LANGUAGE Fortran
    )
    ecbuild_target_fortran_module_directory(
      TARGET ectrans_test_analytic_${platform}_${precision}
      MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests/trans/analytic/module/analytic_test_${platform}_${precision}
    )

    foreach( mpi ${ntasks} )
      set( test_name ectrans_test_analytic_${platform}_mpi${mpi}_${precision} )
      ecbuild_add_test(
        TARGET ${test_name}
        COMMAND ectrans_test_analytic_${platform}_${precision}
        MPI ${mpi}
        OMP 1
      )
    endforeach() # mpi

  endforeach() # platform

endforeach() # precision
