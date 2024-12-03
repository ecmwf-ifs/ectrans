# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

find_package( MPI )

# --------------------------------------------------------------------------------------------------
# Add a test for installation of ecTrans
# --------------------------------------------------------------------------------------------------

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

# --------------------------------------------------------------------------------------------------
# Add a test for SETUP_TRANS0
# --------------------------------------------------------------------------------------------------

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

# --------------------------------------------------------------------------------------------------
# Add a test for tangent-linear/adjoint correspondence (CPU version only)
# --------------------------------------------------------------------------------------------------

if( HAVE_CPU )
  if( HAVE_DOUBLE_PRECISION )
    set( trans trans_dp )
    set( parkind parkind_dp )
  else()
    set( trans trans_sp )
    set( parkind parkind_sp )
  endif()

  ecbuild_add_test(TARGET ectrans_test_adjoint
    SOURCES trans/test_adjoint.F90
    LIBS ${trans} ${parkind}
    LINKER_LANGUAGE Fortran
  )
  if( TEST ectrans_test_adjoint AND HAVE_OMP )
    target_link_libraries( ectrans_test_adjoint OpenMP::OpenMP_Fortran )
  endif()
endif()

# --------------------------------------------------------------------------------------------------
# Add tests for common call patterns of ecTrans, using the benchmark program
# This tests CPU and/or GPU versions, depending on which are enabled
# --------------------------------------------------------------------------------------------------

macro(ectrans_set_test_properties target)
  if( "${target}" MATCHES "gpu" )
    set_tests_properties(${target} PROPERTIES LABELS "gpu;Fortran")
  endif()
endmacro()

# Determine which benchmarks are available
set( benchmarks "" )
if( TARGET ectrans-benchmark-cpu-dp )
  list( APPEND benchmarks ectrans-benchmark-cpu-dp )
endif()
if( TARGET ectrans-benchmark-cpu-sp )
  list( APPEND benchmarks ectrans-benchmark-cpu-sp )
endif()
if( TARGET ectrans-benchmark-gpu-dp )
  list( APPEND benchmarks ectrans-benchmark-gpu-dp )
endif()
if( TARGET ectrans-benchmark-gpu-sp )
  list( APPEND benchmarks ectrans-benchmark-gpu-sp )
endif()

foreach( benchmark ${benchmarks} )
  # Establish which task/thread parameters to test
  set( ntasks 0 )
  set( nthreads 1 )
  if( HAVE_MPI )
    list( APPEND ntasks 1 2 )
  endif()
  if( ${benchmark} MATCHES "cpu" )
    if( HAVE_OMP )
      list( APPEND nthreads 4 8 )
    endif()
  endif()

  # Add test for each parameter combination
  foreach( mpi ${ntasks} )
    foreach( omp ${nthreads} )
      # TCO47 truncation
      set( t 47 )
      set( grid O48 )

      # Base arguments -> 2 iterations, memory consumption/pinning information, spectral norms, and
      # verbose output
      set( base_args "--niter 2 --meminfo --norms -v" )

      set (base_title "${benchmark}_T${t}_${grid}_mpi${mpi}_omp${omp}")

      # Check it works with 0 3D scalar fields
      # This test doesn't work on GPU -> should we delete it?
      if( NOT "${benchmark}" MATCHES "-gpu-" )
        ecbuild_add_test( TARGET ${base_title}_nfld0
            COMMAND ${benchmark}
            ARGS --truncation ${t} --grid ${grid} --nfld 0 --check 100 ${baseargs}
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${base_title}_nfld0 )
      endif()

      # Check it works with 10 3D scalar fields
      ecbuild_add_test( TARGET ${base_title}_nfld10
          COMMAND ${benchmark}
          ARGS --truncation ${t} --grid ${grid} --nfld 10 --check 100 ${baseargs}
          MPI ${mpi}
          OMP ${omp}
      )
      ectrans_set_test_properties( ${base_title}_nfld10 )

      # Check it works with 10 3D scalar fields and 20 levels
      ecbuild_add_test( TARGET ${base_title}_nfld10_nlev20
          COMMAND ${benchmark}
          ARGS --truncation ${t} --grid ${grid} --nfld 10 --nlev 20 --check 100 ${baseargs}
          MPI ${mpi}
          OMP ${omp}
      )
      ectrans_set_test_properties( ${base_title}_nfld10_nlev20 )

      # Check it works with 10 3D scalar fields, 20 levels, and scalar derivatives
      ecbuild_add_test( TARGET ${base_title}_nfld10_nlev20_scders
          COMMAND ${benchmark}
          ARGS --truncation ${t} --grid ${grid} --nfld 10 --nlev 20 --scders --check 100 ${baseargs}
          MPI ${mpi}
          OMP ${omp}
      )
      ectrans_set_test_properties( ${base_title}_nfld10_nlev20_scders )

      # Check it works with 10 3D scalar fields, 20 levels, and wind transforms
      ecbuild_add_test( TARGET ${base_title}_nfld10_nlev20_vordiv
          COMMAND ${benchmark}
          ARGS --truncation ${t} --grid ${grid} --nfld 10 --nlev 20 --vordiv --check 200  ${baseargs}
          MPI ${mpi}
          OMP ${omp}
      )
      ectrans_set_test_properties( ${base_title}_nfld10_nlev20_vordiv )

      # Check it works with 10 3D scalar fields, 20 levels, wind transforms, and wind derivatives
      ecbuild_add_test( TARGET ${base_title}_nfld10_nlev20_vordiv_uvders
          COMMAND ${benchmark}
          ARGS --truncation ${t} --grid ${grid} --nfld 10 --nlev 20 --vordiv --uvders --check 200 ${baseargs}
          MPI ${mpi}
          OMP ${omp}
      )
      ectrans_set_test_properties( ${base_title}_nfld10_nlev20_vordiv_uvders )

      # Check it works with 10 3D scalar fields, 20 levels, and NPROMA=16
      ecbuild_add_test( TARGET ${base_title}_nfld10_nlev20_nproma16
          COMMAND ${benchmark}
          ARGS --truncation ${t} --grid ${grid} --nfld 10 --nlev 20 --nproma 16 --check 100 ${baseargs}
          MPI ${mpi}
          OMP ${omp}
      )
      ectrans_set_test_properties( ${base_title}_nfld10_nlev20_nproma16 )

      if( ${benchmark} MATCHES "cpu" )
        # Check it works with 10 3D scalar fields, 20 levels, and the fast Legendre tranform (CPU only)
        ecbuild_add_test( TARGET ${base_title}_nfld10_nlev20_flt
            COMMAND ${benchmark}
            ARGS --truncation ${t} --grid ${grid} --nfld 10 --nlev 20 --flt --check 4000 ${baseargs}
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${base_title}_nfld10_nlev20_flt )
      endif()
    endforeach()
  endforeach()
endforeach()

# --------------------------------------------------------------------------------------------------
# Add tests for etrans
# --------------------------------------------------------------------------------------------------

if( HAVE_ETRANS )
  foreach( prec dp sp )
    if( TARGET ectrans-lam-benchmark-cpu-${prec} )
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
          set( nlon 48 )
          set( nlat 40 )
          ecbuild_add_test( TARGET ectrans_lam_test_benchmark_${prec}_${nlon}x${nlat}_mpi${mpi}_omp${omp}_nfld0
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS --nlon ${nlon} --nlat ${nlat} --niter 2 --nfld 0 --meminfo --check 100 --norms -v
              MPI ${mpi}
              OMP ${omp}
          )
          ecbuild_add_test( TARGET ectrans_lam_test_benchmark_${prec}_${nlon}x${nlat}_mpi${mpi}_omp${omp}_nfld10
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS --nlon ${nlon} --nlat ${nlat} --niter 2 --nfld 10 --meminfo --check 100 --norms -v
              MPI ${mpi}
              OMP ${omp}
          )
          ecbuild_add_test( TARGET ectrans_lam_test_benchmark_${prec}_${nlon}x${nlat}_mpi${mpi}_omp${omp}_nfld10_nlev20
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS --nlon ${nlon} --nlat ${nlat} --niter 2 --nfld 10 --nlev 20 --check 100 --norms -v
              MPI ${mpi}
              OMP ${omp}
          )
          ecbuild_add_test( TARGET ectrans_lam_test_benchmark_${prec}_${nlon}x${nlat}_mpi${mpi}_omp${omp}_nfld10_nlev20_scders
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS --nlon ${nlon} --nlat ${nlat} --niter 2 --nfld 10 --nlev 20 --scders --check 100 --norms -v
              MPI ${mpi}
              OMP ${omp}
          )
          ecbuild_add_test( TARGET ectrans_lam_test_benchmark_${prec}_${nlon}x${nlat}_mpi${mpi}_omp${omp}_nfld10_nlev20_vordiv
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS --nlon ${nlon} --nlat ${nlat} --niter 2 --nfld 10 --nlev 20 --vordiv --check 100 --norms -v
              MPI ${mpi}
              OMP ${omp}
          )
          ecbuild_add_test( TARGET ectrans_lam_test_benchmark_${prec}_${nlon}x${nlat}_mpi${mpi}_omp${omp}_nfld10_nlev20_vordiv_uvders
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS --nlon ${nlon} --nlat ${nlat} --niter 2 --nfld 10 --nlev 20 --vordiv --uvders --check 100 --norms -v
              MPI ${mpi}
              OMP ${omp}
          )
          ecbuild_add_test( TARGET ectrans_lam_test_benchmark_${prec}_${nlon}x${nlat}_mpi${mpi}_omp${omp}_nfld10_nlev20_nproma16
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS --nlon ${nlon} --nlat ${nlat} --niter 2 --nfld 10 --nlev 20 --nproma 16 --check 100 --norms -v
              MPI ${mpi}
              OMP ${omp}
          )
        endforeach()
      endforeach()
    endif()
  endforeach()

endif()

# --------------------------------------------------------------------------------------------------
# Add tests for transi
# --------------------------------------------------------------------------------------------------

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

  if( HAVE_GPU )
    ecbuild_add_library( TARGET ectrans_test_gpu
      SOURCES     transi/transi_test.h transi/transi_test.c
      PUBLIC_LIBS transi_gpu_dp
      NOINSTALL
    )
    target_compile_definitions( ectrans_test PUBLIC ${TEST_DEFINITIONS} )
  endif()

  ecbuild_add_test( TARGET ectrans_test_transi_program
    SOURCES   transi/transi_test_program.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    ENVIRONMENT TRANS_USE_MPI=0 )

  if( HAVE_GPU )
    ecbuild_add_test( TARGET ectrans_test_transi_program_gpu
      SOURCES   transi/transi_test_program.c
      LIBS      ectrans_test_gpu
      LINKER_LANGUAGE C
      DEFINITIONS GPU_VERSION
      ENVIRONMENT TRANS_USE_MPI=0 )
    set_tests_properties(ectrans_test_transi_program_gpu PROPERTIES LABELS gpu)
  endif()

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

  ecbuild_add_option( FEATURE MEMORY_TESTS DEFAULT ON DESCRIPTION "Enable memory tests" )
  if( NOT HAVE_MEMORY_TESTS )
    set_tests_properties( ectrans_test_transi_memory ectrans_test_transi_memory_lonlat PROPERTIES DISABLED ON )
  endif()
endif()
