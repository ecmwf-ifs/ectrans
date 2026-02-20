# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Tests for adjoint functionality

# --------------------------------------------------------------------------------------------------
# Add a test for tangent-linear/adjoint correspondence of INV_TRANS and DIR_TRANS combined (CPU
# version only)
# --------------------------------------------------------------------------------------------------

if( HAVE_CPU )
  foreach( precision ${precisions} )
    foreach( mpi ${ntasks} )
      set( test_name ectrans_test_adjoint_mpi${mpi}_${precision} )
      ecbuild_add_test(TARGET ${test_name}
        SOURCES test_adjoint.F90
        LIBS trans_${precision} parkind_${precision}
        LINKER_LANGUAGE Fortran
        MPI ${mpi}
      )
      if( TEST ${test_name} AND HAVE_OMP )
        target_link_libraries( ${test_name} OpenMP::OpenMP_Fortran )
      endif()
    endforeach()
  endforeach()
endif()

# --------------------------------------------------------------------------------------------------
# Add tests for tangent-linear/adjoint correspondence of INV_TRANS and DIR_TRANS
# --------------------------------------------------------------------------------------------------

foreach( platform ${platforms} )
  
  if( "${platform}" MATCHES "cpu" )
    set( platform_tag "" )
  endif()
  if( "${platform}" MATCHES "gpu" )
    set( platform_tag "_gpu" )
  endif()
  
  if( HAVE_DOUBLE_PRECISION )
    set( trans trans${platform_tag}_dp )
    set( parkind parkind_dp )
  else()
    set( trans trans${platform_tag}_sp )
    set( parkind parkind_sp )
  endif()

  foreach( mpi ${ntasks} )
    ecbuild_add_test(TARGET ectrans_test_dirtrans_adjoint_${platform}_mpi${mpi}
      SOURCES test_dirtrans_adjoint.F90
      LIBS ${trans} ${parkind}
      LINKER_LANGUAGE Fortran
      MPI ${mpi}
      OMP 1
    )
    if( TEST ectrans_test_dirtrans_adjoint_${platform}_mpi${mpi} AND HAVE_OMP )
      target_link_libraries( ectrans_test_dirtrans_adjoint_${platform}_mpi${mpi} OpenMP::OpenMP_Fortran )
    endif()
    set_tests_properties(ectrans_test_dirtrans_adjoint_${platform}_mpi${mpi} PROPERTIES LABELS "${platform};Fortran")

    ecbuild_add_test(TARGET ectrans_test_invtrans_adjoint_${platform}_mpi${mpi}
      SOURCES test_invtrans_adjoint.F90
      LIBS ${trans} ${parkind}
      LINKER_LANGUAGE Fortran
      MPI ${mpi}
      OMP 1
    )
    if( TEST ectrans_test_invtrans_adjoint_${platform}_mpi${mpi} AND HAVE_OMP )
      target_link_libraries( ectrans_test_invtrans_adjoint_${platform}_mpi${mpi} OpenMP::OpenMP_Fortran )
    endif()
    set_tests_properties(ectrans_test_invtrans_adjoint_${platform}_mpi${mpi} PROPERTIES LABELS "${platform};Fortran")
    
  endforeach() # mpi comm size

endforeach() # platform

# --------------------------------------------------------------------------------------------------
# Add a test for tangent-linear/adjoint correspondence of GPNORM_TRANSTL/AD (CPU version only)
# --------------------------------------------------------------------------------------------------

if( HAVE_CPU )
  foreach( precision ${precisions} )
    foreach( mpi ${ntasks} )
      set( test_name ectrans_test_gpnorm_trans_adjoint_mpi${mpi}_${precision} )
      ecbuild_add_test(TARGET ${test_name}
        SOURCES test_gpnorm_adjoint.F90
        LIBS trans_${precision} parkind_${precision}
        LINKER_LANGUAGE Fortran
        MPI ${mpi}
      )
      if( TEST ${test_name} AND HAVE_OMP )
        target_link_libraries( ${test_name} OpenMP::OpenMP_Fortran )
      endif()
    endforeach()
  endforeach()
endif()
