# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### trans_test_install

if( HAVE_TESTS )

  find_package( MPI )
  configure_file( test-install.sh.in ${CMAKE_CURRENT_BINARY_DIR}/test-install.sh @ONLY )

  unset( _test_args )
  if( CMAKE_TOOLCHAIN_FILE )
    list( APPEND _test_args "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}" )
  endif()
  foreach( lang C CXX Fortran )
    if( CMAKE_${lang}_COMPILER )
      list( APPEND _test_args "-DCMAKE_${lang}_COMPILER=${CMAKE_${lang}_COMPILER}" )
    endif()
  endforeach()
  if( NOT HAVE_DOUBLE_PRECISION )
    list( APPEND _test_args "-DCOMPONENTS=single" )
  endif() 

  add_test( NAME ectrans_test_install
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test-install.sh ${_test_args} )
endif()


if( HAVE_DOUBLE_PRECISION )
  set( trans trans_dp )
else()
  set( trans trans_sp )
endif()

ecbuild_add_test(TARGET ectrans_test_adjoint
  SOURCES trans/test_adjoint.F90
  LIBS ${trans}
  LINKER_LANGUAGE Fortran
)
if( TEST ectrans_test_adjoint AND HAVE_OMP )
  target_link_libraries( ectrans_test_adjoint OpenMP::OpenMP_Fortran )
endif()

set( ntasks 0 )
set( nthreads 1 )
if( HAVE_MPI )
  list( APPEND ntasks 1 2 )
endif()
if( HAVE_OMP )
  list( APPEND nthreads 4 8 )
endif()

foreach( prec dp sp )
  foreach( mpi ${ntasks} )
    foreach( omp ${nthreads} )
        set( t 47 )
        set( grid O48 )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld0
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 0 -v --meminfo --check 100
            MPI ${mpi}
            OMP ${omp}
        )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --meminfo --check 100
            MPI ${mpi}
            OMP ${omp}
        )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --check 100
            MPI ${mpi}
            OMP ${omp}
        )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_scders
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --scders --check 100
            MPI ${mpi}
            OMP ${omp}
        )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_vordiv
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --vordiv --check 100
            MPI ${mpi}
            OMP ${omp}
        )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_vordiv_uvders
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --vordiv --uvders --check 100
            MPI ${mpi}
            OMP ${omp}
        )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_flt
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --flt --check 1000
            MPI ${mpi}
            OMP ${omp}
        )
        ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_nproma16
            COMMAND ectrans-benchmark-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --nproma 16 --check 100
            MPI ${mpi}
            OMP ${omp}
        )
    endforeach()
  endforeach()
endforeach()


if( HAVE_TRANSI )

  check_include_files( malloc.h       EC_HAVE_MALLOC_H      )
  ecbuild_debug_var( EC_HAVE_MALLOC_H )


  if( EC_HAVE_MALLOC_H )
    list( APPEND  TEST_DEFINITIONS
      TRANSI_HAVE_MEMORY
    )
  else()
    ecbuild_warn( "ectrans tests checking memory leaks are disabled as malloc.h was not found" )
  endif()

  ecbuild_add_library( TARGET ectrans_test
    SOURCES     transi/transi_test.h transi/transi_test.c
    PUBLIC_LIBS transi_dp
    NOINSTALL
  )
  target_compile_definitions( ectrans_test PUBLIC ${TEST_DEFINITIONS} )

  ecbuild_add_test( TARGET ectrans_test_transi_program
    SOURCES   transi/transi_test_program.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_timings
    SOURCES   transi/transi_test_timings.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_lonlat
    SOURCES   transi/transi_test_lonlat.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_io
    SOURCES   transi/transi_test_io.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_memory
    SOURCES   transi/transi_test_memory.c
    LIBS      ectrans_test
    CONDITION EC_HAVE_MALLOC_H
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_memory_lonlat
    SOURCES   transi/transi_test_memory_lonlat.c
    LIBS      ectrans_test
    CONDITION EC_HAVE_MALLOC_H
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_vordiv_to_UV
    SOURCES   transi/transi_test_vordiv_to_UV.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_dirtrans_adjoint
    SOURCES   transi/transi_test_dirtrans_adjoint.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_invtrans_adjoint
    SOURCES   transi/transi_test_invtrans_adjoint.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_lonlat_diff_incr
    SOURCES   transi/transi_test_lonlat_diff_incr.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  if( HAVE_TESTS )
    ecbuild_add_option( FEATURE MEMORY_TESTS DEFAULT ON DESCRIPTION "Enable memory tests" )
    if( NOT HAVE_MEMORY_TESTS )
      set_tests_properties( ectrans_test_transi_memory ectrans_test_transi_memory_lonlat PROPERTIES DISABLED ON )
    endif()
  endif()

endif()
