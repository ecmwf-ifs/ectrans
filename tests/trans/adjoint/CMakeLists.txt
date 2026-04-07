# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Tests for adjoint functionality

foreach( precision ${precisions} )

  ecbuild_add_library(
    TARGET           adjoint_utils_${precision}
    LINKER_LANGUAGE  Fortran
    SOURCES          utils.F90
    PUBLIC_LIBS      fiat trans_${precision}
  )
  ecbuild_target_fortran_module_directory(
    TARGET adjoint_utils_${precision}
    MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests/trans/adjoint/module/adjoint_utils_${precision}
  )

  foreach( platform ${platforms} )

    if( "${platform}" MATCHES "cpu" )
      set( platform_tag "" )
    endif()
    if( "${platform}" MATCHES "gpu" )
      set( platform_tag "_gpu" )
    endif()

    foreach( mpi ${ntasks} )
      # Add a test for tangent-linear/adjoint correspondence of INV_TRANS and DIR_TRANS combined
      ecbuild_add_test(TARGET ectrans_test_adjoint_${platform}_mpi${mpi}_${precision}
        SOURCES test_adjoint.F90
        LIBS trans${platform_tag}_${precision} parkind_${precision} adjoint_utils_${precision}
        LINKER_LANGUAGE Fortran
        MPI ${mpi}
        OMP 1
      )

      # Add test for tangent-linear/adjoint correspondence of DIR_TRANS only
      ecbuild_add_test(TARGET ectrans_test_dirtrans_adjoint_${platform}_mpi${mpi}_${precision}
        SOURCES test_dirtrans_adjoint.F90
        LIBS trans${platform_tag}_${precision} parkind_${precision} adjoint_utils_${precision}
        LINKER_LANGUAGE Fortran
        MPI ${mpi}
        OMP 1
      )

      # Add test for tangent-linear/adjoint correspondence of INV_TRANS only
      ecbuild_add_test(TARGET ectrans_test_invtrans_adjoint_${platform}_mpi${mpi}_${precision}
        SOURCES test_invtrans_adjoint.F90
        LIBS trans${platform_tag}_${precision} parkind_${precision} adjoint_utils_${precision}
        LINKER_LANGUAGE Fortran
        MPI ${mpi}
        OMP 1
      )

      # GPNORMTRANS_TL/AD are not yet implemented for GPU
      if( NOT "${platform}" MATCHES "gpu" )
        # Add a test for tangent-linear/adjoint correspondence of GPNORM_TRANSTL/AD
        ecbuild_add_test(TARGET ectrans_test_gpnorm_trans_adjoint_${platform}_mpi${mpi}_${precision}
          SOURCES test_gpnorm_adjoint.F90
          LIBS trans${platform_tag}_${precision} parkind_${precision} adjoint_utils_${precision}
          LINKER_LANGUAGE Fortran
          MPI ${mpi}
          OMP 1
        )
      endif()

    endforeach() # mpi

  endforeach() # platform

endforeach() # precision
