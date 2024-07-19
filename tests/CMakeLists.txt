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
  foreach( lang C CXX Fortran )
    if( CMAKE_${lang}_FLAGS )
      list( APPEND _test_args "-DCMAKE_${lang}_FLAGS=${CMAKE_${lang}_FLAGS}" )
    endif()
  endforeach()
  if( CMAKE_EXE_LINKER_FLAGS )
    list( APPEND _test_args "-DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}" )
  endif()
  if( NOT HAVE_DOUBLE_PRECISION )
    list( APPEND _test_args "-DCOMPONENTS=single" )
  endif()

  add_test( NAME ectrans_test_install
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test-install.sh ${_test_args} )

  if( HAVE_DOUBLE_PRECISION )
    set( trans trans_dp )
    set( parkind parkind_dp )
  else()
    set( trans trans_sp )
    set( parkind parkind_sp )
  endif()

  ecbuild_add_executable(
    TARGET  ectrans_test_setup_trans0
    SOURCES trans/test_setup_trans0.F90
    LIBS    ectrans_common
    LINKER_LANGUAGE Fortran
    NOINSTALL)
  set( ntasks 0 )
  if( HAVE_MPI )
    list( APPEND ntasks 1 2 )
  endif()
  foreach( mpi ${ntasks} )
      ecbuild_add_test( TARGET ectrans_test_setup_trans0_mpi${mpi}
          COMMAND ectrans_test_setup_trans0
          MPI ${mpi}
      )
  endforeach()

  ecbuild_add_test(TARGET ectrans_test_adjoint
    SOURCES trans/test_adjoint.F90
    LIBS ${trans} ${parkind}
    LINKER_LANGUAGE Fortran
  )
  if( TEST ectrans_test_adjoint AND HAVE_OMP )
    target_link_libraries( ectrans_test_adjoint OpenMP::OpenMP_Fortran )
  endif()


  foreach( prec dp sp )
    if( TARGET ectrans-benchmark-cpu-${prec} )
      set( ntasks 0 )
      set( nthreads 1 )
      if( HAVE_MPI )
        list( APPEND ntasks 1 2 )
      endif()
      if( HAVE_OMP )
        list( APPEND nthreads 4 8 )
      endif()
      foreach( mpi ${ntasks} )
        foreach( omp ${nthreads} )
          foreach( callmode 1 2 )
            set( t 47 )
            set( grid O48 )
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld0_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 0 --meminfo --check 100 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --meminfo --check 100 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --check 100 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_scders_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --scders --check 100 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_vordiv_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --vordiv --check 100 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_vordiv_uvders_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --vordiv --uvders --check 100 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
            # TODO: Find out why the FLT gives so much higher errors
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_flt_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --flt --check 1000000 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
            ecbuild_add_test( TARGET ectrans_test_benchmark_${prec}_T${t}_${grid}_mpi${mpi}_omp${omp}_nfld10_nlev20_nproma16_callmode${callmode}
                COMMAND ectrans-benchmark-cpu-${prec} ARGS --truncation ${t} --grid ${grid} --niter 2 --nfld 10 --nlev 20 --nproma 16 --check 100 --norms -v --callmode ${callmode}
                MPI ${mpi}
                OMP ${omp}
            )
          endforeach()
        endforeach()
      endforeach()
    endif()
  endforeach()

endif()

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
