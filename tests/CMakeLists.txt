# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# --------------------------------------------------------------------------------------------------
# First establish which precisions and platforms to test
# --------------------------------------------------------------------------------------------------

set( precisions "" )
if( HAVE_DOUBLE_PRECISION )
  list( APPEND precisions "dp" )
endif()
if( HAVE_SINGLE_PRECISION )
  list( APPEND precisions "sp" )
endif()

set( platforms )
if( HAVE_CPU )
  list( APPEND platforms cpu )
endif()
if( HAVE_GPU )
  list( APPEND platforms gpu )
endif()

# --------------------------------------------------------------------------------------------------
# Then set the MPI and OpenMP testing configurations
# --------------------------------------------------------------------------------------------------

set( ntasks 0 )
if( HAVE_MPI )
  list( APPEND ntasks 1 4 )
endif()
set( nthreads 1 )
if( HAVE_OMP )
  list( APPEND nthreads 8 )
endif()

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
# Add API tests
# --------------------------------------------------------------------------------------------------

add_subdirectory( trans/api )

# --------------------------------------------------------------------------------------------------
# Add adjoint tests
# --------------------------------------------------------------------------------------------------

add_subdirectory( trans/adjoint )

# --------------------------------------------------------------------------------------------------
# Add field-view metadata test
# --------------------------------------------------------------------------------------------------

if( HAVE_CPU )
  set( field_view_platform_tag "" )
else()
  set( field_view_platform_tag "gpu_" )
endif()
if( HAVE_DOUBLE_PRECISION )
  set( field_view_precision dp )
else()
  set( field_view_precision sp )
endif()

ecbuild_add_test( TARGET ectrans_test_field_view
  SOURCES trans/test_field_view.F90
  LIBS ectrans_field_view_${field_view_platform_tag}${field_view_precision}
  LINKER_LANGUAGE Fortran )

# --------------------------------------------------------------------------------------------------
# Add tests for common call patterns of ecTrans, using the benchmark program
# This tests CPU and/or GPU versions, depending on which are enabled
# --------------------------------------------------------------------------------------------------

set( checksums_dir ${CMAKE_CURRENT_BINARY_DIR}/checksums )
file( MAKE_DIRECTORY ${checksums_dir} )
set( local_reference_dir ${CMAKE_CURRENT_BINARY_DIR}/local-reference )
file( MAKE_DIRECTORY ${local_reference_dir} )
configure_file( ectrans-update-references.sh.in
                ${CMAKE_CURRENT_BINARY_DIR}/ectrans-update-references.sh @ONLY )

set( ectrans_update_reference_benchmarks )

set( reference_dir ${CMAKE_CURRENT_BINARY_DIR}/reference )
if( IS_SYMLINK "${reference_dir}" )
  file( REMOVE "${reference_dir}" )
elseif( EXISTS "${reference_dir}" )
  ecbuild_critical( "Cannot create reference symlink because path already exists: ${reference_dir}" )
endif()
file( CREATE_LINK "${ECTRANS_TEST_REFERENCE_DIR}" "${reference_dir}" SYMBOLIC )

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
  if( benchmark MATCHES "ectrans-benchmark-cpu" )
    list( APPEND ectrans_update_reference_benchmarks ${benchmark} )
    set( source_prefix "${checksums_dir}/${benchmark}_T47_O48_mpi0_omp1_callmode1" )

    ectrans_checksum_reference_name( reference_name ${benchmark} T47 O48 nfld0 nlev1 niter2 scders0 uvders0 vordiv0 flt0 )
    ectrans_add_local_reference_file( ${source_prefix}_nfld0 ${reference_name} )

    ectrans_checksum_reference_name( reference_name ${benchmark} T47 O48 nfld10 nlev20 niter2 scders0 uvders0 vordiv0 flt0 )
    ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20 ${reference_name} )

    ectrans_checksum_reference_name( reference_name ${benchmark} T47 O48 nfld10 nlev20 niter2 scders1 uvders1 vordiv0 flt0 )
    ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20_derivatives ${reference_name} )

    ectrans_checksum_reference_name( reference_name ${benchmark} T47 O48 nfld10 nlev20 niter2 scders0 uvders0 vordiv1 flt0 )
    ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20_vordiv ${reference_name} )

    ectrans_checksum_reference_name( reference_name ${benchmark} T47 O48 nfld10 nlev20 niter2 scders0 uvders0 vordiv0 flt1 )
    ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20_flt ${reference_name} )

  endif()

  # Establish which task/thread parameters to test
  set( ntasks 0 )
  set( nthreads 1 )
  if( HAVE_MPI )
    list( APPEND ntasks 1 4 )
  endif()
  if( ${benchmark} MATCHES "cpu" )
    if( HAVE_OMP )
      list( APPEND nthreads 8 )
    endif()
  endif()

  # Add test for each parameter combination
  foreach( mpi ${ntasks} )
    foreach( omp ${nthreads} )
      # This test suite assumes a 4-core machine (e.g. Github runner) so we don't want to
      # oversubscribe the cores
      if( "${mpi}" STREQUAL "4" AND ${omp} GREATER "1" )
        continue()
      endif()

      foreach( callmode 1 2 )

        # Base arguments -> 2 iterations, memory consumption/pinning information, spectral norms, and
        # verbose output
        set( t 47 )
        set( grid O48 )
        set( niter 2 )
        set( base_args --truncation ${t} --grid ${grid} --callmode ${callmode} --niter ${niter} --niter-checksums ${niter} --meminfo --norms -v )
        set (base_title "${benchmark}_T${t}_${grid}_mpi${mpi}_omp${omp}_callmode${callmode}")

        # Check it works with 0 3D scalar fields
        set( test_name ${base_title}_nfld0 )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 0 --check 100 --callmode ${callmode} --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld0 nlev1 niter${niter} scders0 uvders0 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        # Check it works with 10 3D scalar fields and 20 levels
        set( test_name ${base_title}_nfld10_nlev20 )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 10 --nlev 20 --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        # Check it works with 10 3D scalar fields, 20 levels, and derivatives
        set( test_name ${base_title}_nfld10_nlev20_derivatives )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 10 --nlev 20 --scders --uvders --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders1 uvders1 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )
        # Check it works with 10 3D scalar fields, 20 levels, and vordiv in grid point space
        set( test_name ${base_title}_nfld10_nlev20_vordiv )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 10 --nlev 20 --vordiv --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv1 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        # Check it works with 10 3D scalar fields, 20 levels, and NPROMA=16
        set( test_name ${base_title}_nfld10_nlev20_nproma16 )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 10 --nlev 20 --nproma 16 --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        # Check it works with LALLOPERM=F (cpu only)
        if( ${benchmark} MATCHES "cpu" )
          set( test_name ${base_title}_lalloperm_f )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ${benchmark}
              ARGS ${base_args} --nfld 10 --nlev 20 --deallocate-foubuf-temps --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )
        endif()

        # Check it works with 10 3D scalar fields, 20 levels, and NPROMATR=20 (feature only works for CPU)
        if( ${benchmark} MATCHES "cpu" AND ${callmode} MATCHES "1" )
          set( test_name ${base_title}_nfld10_nlev20_npromatr20 )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ${benchmark}
              ARGS ${base_args} --nfld 10 --nlev 20 --npromatr 20 --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )
        endif()

        if( ${benchmark} MATCHES "cpu" )
          # Check it works with 10 3D scalar fields, 20 levels, and the fast Legendre tranform (CPU only)
          # TODO: Find out why the FLT gives so much higher errors
          set( test_name ${base_title}_nfld10_nlev20_flt )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ${benchmark}
              ARGS ${base_args} --nfld 10 --nlev 20 --flt --check 1000000 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt1 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )
        endif()
      endforeach()
    endforeach()
  endforeach()
endforeach()

# --------------------------------------------------------------------------------------------------
# Add tests for field API interface, using the benchmark program
# --------------------------------------------------------------------------------------------------

if (HAVE_FIELD_API)
foreach( benchmark ${benchmarks} )
  # Establish which task/thread parameters to test
  set( ntasks 0 )
  set( nthreads 1 )
  if( HAVE_MPI )
    list( APPEND ntasks 2 )
  endif()
  if( ${benchmark} MATCHES "cpu" )
    if( HAVE_OMP )
      list( APPEND nthreads 8)
    endif()
  endif()

  # Add test for each parameter combination
  foreach( mpi ${ntasks} )
    foreach( omp ${nthreads} )
      foreach( callmode 1 2 )

        set( t 47 )
        set( grid O48 )
        set( niter 2 )

        set( base_args --truncation ${t} --grid ${grid} --niter ${niter} --niter-checksums ${niter} --callmode ${callmode} --meminfo --norms -v --field-api)
        set (base_title "${benchmark}_T${t}_${grid}_mpi${mpi}_omp${omp}_callmode${callmode}_field_api")

        #----
        # Check it works with 0 3D scalar fields
        set( test_name ${base_title}_nfld0 )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 0 --check 100 --callmode ${callmode} --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld0 nlev1 niter${niter} scders0 uvders0 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )
        #----
        # Check it works with 10 3D scalar fields and 20 levels
        set( test_name ${base_title}_nfld10_nlev20 )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 10 --nlev 20 --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        #----
        # Check it works with 10 3D scalar fields, 20 levels, and derivatives
        set( test_name ${base_title}_nfld10_nlev20_derivatives )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
          ARGS ${base_args} --nfld 10 --nlev 20 --scders --uvders --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders1 uvders1 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        #----
        # Check it works with 10 3D scalar fields, 20 levels, and wind transforms
        set( test_name ${base_title}_nfld10_nlev20_vordiv )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
            ARGS ${base_args} --nfld 10 --nlev 20 --vordiv --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
            MPI ${mpi}
            OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv1 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        #----
        # Check it works with 10 3D scalar fields, 20 levels, and NPROMA=16
        set( test_name ${base_title}_nfld10_nlev20_nproma16 )
        ecbuild_add_test( TARGET ${test_name}
            COMMAND ${benchmark}
          ARGS ${base_args} --nfld 10 --nlev 20 --nproma 16 --check 100 --dump-checksums ${checksums_dir}/${test_name}.checksums
          MPI ${mpi}
          OMP ${omp}
        )
        ectrans_set_test_properties( ${test_name} )
        ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
        ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        #----
        if( ${benchmark} MATCHES "cpu" )
          # Check it works with 10 3D scalar fields, 20 levels, and the fast Legendre tranform (CPU only)
          set( test_name ${base_title}_nfld10_nlev20_flt )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ${benchmark}
              ARGS ${base_args} --nfld 10 --nlev 20 --flt --check 1000000 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} T${t} ${grid} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt1 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )
        endif()

      endforeach() # callmode
    endforeach() # omp
  endforeach() # mpi
endforeach() # benchmark
endif(HAVE_FIELD_API)

# --------------------------------------------------------------------------------------------------
# Add tests for common call patterns of ecTrans LAM benchmark (i.e. etrans), using the benchmark
# program
# --------------------------------------------------------------------------------------------------

if( HAVE_ETRANS )
  # Set resolution
  set( nlon 48 )
  set( nlat 40 )
  set( niter 2 )

  foreach( prec dp sp )
    if( TARGET ectrans-lam-benchmark-cpu-${prec} )
      set( benchmark ectrans-lam-benchmark-cpu-${prec} )
      list( APPEND ectrans_update_reference_benchmarks ${benchmark} )
      set( source_prefix "${checksums_dir}/${benchmark}_${nlon}x${nlat}_mpi0_omp1" )

      ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld0 nlev1 niter${niter} scders0 uvders0 vordiv0 flt0 )
      ectrans_add_local_reference_file( ${source_prefix}_nfld0 ${reference_name} )

      ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev1 niter${niter} scders0 uvders0 vordiv0 flt0 )
      ectrans_add_local_reference_file( ${source_prefix}_nfld10 ${reference_name} )

      ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
      ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20 ${reference_name} )

      ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders1 uvders0 vordiv0 flt0 )
      ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20_scders ${reference_name} )

      ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv1 flt0 )
      ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20_vordiv ${reference_name} )

      ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders1 vordiv1 flt0 )
      ectrans_add_local_reference_file( ${source_prefix}_nfld10_nlev20_vordiv_uvders ${reference_name} )

      set( nthreads 1 )
      if( HAVE_OMP )
        list( APPEND nthreads 8 )
      endif()

      # Base arguments -> nlat x nlon, 2 iterations, memory consumption/pinning information,
      # spectral norms, and verbose output
      set( base_args --nlon ${nlon} --nlat ${nlat} --niter ${niter} --niter-checksums ${niter} --meminfo --norms -v )

      foreach( mpi ${ntasks} )
        foreach( omp ${nthreads} )
          # This test suite assumes a 4-core machine (e.g. Github runner) so we don't want to
          # oversubscribe the cores
          if( "${mpi}" STREQUAL "4" AND "${omp}" GREATER "1" )
            continue()
          endif()

          set(base_title "${benchmark}_${nlon}x${nlat}_mpi${mpi}_omp${omp}" )

          set( test_name ${base_title}_nfld0 )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 0 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld0 nlev1 niter${niter} scders0 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

          set( test_name ${base_title}_nfld10 )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 10 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev1 niter${niter} scders0 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

          set( test_name ${base_title}_nfld10_nlev20 )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 10 --nlev 20 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

          set( test_name ${base_title}_nfld10_nlev20_scders )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 10 --nlev 20 --scders --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties(  ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders1 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

          set( test_name ${base_title}_nfld10_nlev20_vordiv )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 10 --nlev 20 --vordiv --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv1 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

          set( test_name ${base_title}_nfld10_nlev20_vordiv_uvders )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 10 --nlev 20 --vordiv --uvders --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders1 vordiv1 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

          set( test_name ${base_title}_nfld10_nlev20_nproma16 )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 10 --nlev 20 --nproma 16 --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

          # Check it works with LALLOPERM=F (cpu only)
          set( test_name ${base_title}_lalloperm_f )
          ecbuild_add_test( TARGET ${test_name}
              COMMAND ectrans-lam-benchmark-cpu-${prec} ARGS ${base_args} --nfld 10 --nlev 20 --nproma 16 --deallocate-foubuf-temps --dump-checksums ${checksums_dir}/${test_name}.checksums
              MPI ${mpi}
              OMP ${omp}
          )
          ectrans_set_test_properties( ${test_name} )
          ectrans_checksum_reference_name( reference_name ${benchmark} LAM ${nlon}x${nlat} nfld10 nlev20 niter${niter} scders0 uvders0 vordiv0 flt0 )
          ectrans_add_checksum_reference_tests( ${test_name} ${reference_name} )

        endforeach()
      endforeach()
    endif()
  endforeach()

endif()

get_property( ectrans_update_reference_tests GLOBAL PROPERTY ECTRANS_LOCAL_REFERENCE_SOURCE_TARGETS )
if( ectrans_update_reference_tests )
  list( REMOVE_DUPLICATES ectrans_update_reference_tests )
  list( REMOVE_DUPLICATES ectrans_update_reference_benchmarks )
  string( REPLACE ";" "$|^" ectrans_update_reference_regex "${ectrans_update_reference_tests}" )
  set( ectrans_update_reference_regex "^${ectrans_update_reference_regex}$" )

  add_custom_target( ectrans-update-references
    COMMAND "${CMAKE_CTEST_COMMAND}" --test-dir "${CMAKE_BINARY_DIR}" --output-on-failure --tests-regex "${ectrans_update_reference_regex}"
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/ectrans-update-references.sh"
    DEPENDS ${ectrans_update_reference_benchmarks}
    USES_TERMINAL
    VERBATIM )
  set_target_properties( ectrans-update-references PROPERTIES EXCLUDE_FROM_ALL TRUE )
endif()

# --------------------------------------------------------------------------------------------------
# Add test for split-communicator compatibility of ecTrans
# --------------------------------------------------------------------------------------------------

# It doesn't really matter which precision we use for this test, as long as we use one that's
# actually available
if( HAVE_SINGLE_PRECISION )
  set( trans_lib trans_sp )
  set( parkind_lib parkind_sp )
else()
  set( trans_lib trans_dp )
  set( parkind_lib parkind_dp )
endif()
if( HAVE_CPU AND HAVE_MPI )
  if (fiat_VERSION VERSION_LESS 2.0 )
    ecbuild_warn("Skipping test ectrans_test_split_mpi_comm for fiat version < 2.0; see https://github.com/ecmwf-ifs/fiat/pull/72")
  else()
    ecbuild_add_test(TARGET ectrans_test_split_mpi_comm
      SOURCES trans/test_split_mpi_comm.F90
      LIBS ${trans_lib} ${parkind_lib}
      LINKER_LANGUAGE Fortran
      MPI 4
      OMP 1
    )
  endif()

  if (fiat_VERSION VERSION_LESS 1.6.3 OR (fiat_VERSION VERSION_EQUAL 2.0.0 AND NOT HAVE_FIAT_DEVELOP))
    # Once fiat develop branch is released in version > 2.0.0, we can drop "AND NOT HAVE FIAT_DEVELOP"
    ecbuild_warn("Skipping test ectrans_test_init_mpi_ranks for fiat version (${fiat_VERSION}) < 1.6.3 or "
                 " version == 2.0.0; see https://github.com/ecmwf-ifs/ectrans/pull/387")
  else()
    ecbuild_add_test(TARGET ectrans_test_init_mpi_ranks
      SOURCES trans/test_init_mpi_ranks.F90
      LIBS ${trans_lib} ${parkind_lib}
      LINKER_LANGUAGE Fortran
      MPI 4
      OMP 1
    )
  endif()
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

  find_package( CMath )
  ecbuild_add_library( TARGET ectrans_test
    SOURCES     transi/transi_test.h transi/transi_test.c
    PUBLIC_LIBS transi_dp ${CMATH_LIBRARIES}
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

    ecbuild_add_test( TARGET ectrans_test_transi_invtrans_adjoint_gpu
      SOURCES   transi/transi_test_invtrans_adjoint.c
      LIBS      ectrans_test_gpu
      LINKER_LANGUAGE C
      DEFINITIONS GPU_VERSION
      ENVIRONMENT TRANS_USE_MPI=0 )
    set_tests_properties(ectrans_test_transi_invtrans_adjoint_gpu PROPERTIES LABELS gpu)

    ecbuild_add_test( TARGET ectrans_test_transi_dirtrans_adjoint_gpu
      SOURCES   transi/transi_test_dirtrans_adjoint.c
      LIBS      ectrans_test_gpu
      LINKER_LANGUAGE C
      DEFINITIONS GPU_VERSION
      ENVIRONMENT TRANS_USE_MPI=0 )
    set_tests_properties(ectrans_test_transi_dirtrans_adjoint_gpu PROPERTIES LABELS gpu)
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

  if( NOT HAVE_MEMORY_TESTS )
    set_tests_properties( ectrans_test_transi_memory ectrans_test_transi_memory_lonlat PROPERTIES DISABLED ON )
  endif()

  ecbuild_add_test( TARGET ectrans_test_transi_lam
    SOURCES   transi/transi_test_lam.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    CONDITION HAVE_ETRANS
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_lam_dirtrans_adjoint
    SOURCES   transi/transi_test_lam_dirtrans_adjoint.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    CONDITION HAVE_ETRANS
    ENVIRONMENT TRANS_USE_MPI=0 )

  ecbuild_add_test( TARGET ectrans_test_transi_lam_invtrans_adjoint
    SOURCES   transi/transi_test_lam_invtrans_adjoint.c
    LIBS      ectrans_test
    LINKER_LANGUAGE C
    CONDITION HAVE_ETRANS
    ENVIRONMENT TRANS_USE_MPI=0 )

  if( HAVE_MPI )
    # Test ectrans_test_transi_split_comm requires MPI C
    find_package( MPI COMPONENTS C )
    if ( MPI_C_FOUND )
      ecbuild_add_test( TARGET ectrans_test_transi_split_comm
        SOURCES   transi/transi_test_split_comm.c
        LIBS      ectrans_test MPI::MPI_C
        MPI 2
        LINKER_LANGUAGE C
        ENVIRONMENT TRANS_USE_MPI=1 )
    endif()
  endif()

  if( HAVE_ETRANS_GPU )
    ecbuild_add_test( TARGET ectrans_test_transi_lam_gpu
      SOURCES   transi/transi_test_lam.c
      LIBS      ectrans_test_gpu
      LINKER_LANGUAGE C
      DEFINITIONS GPU_VERSION
      CONDITION HAVE_ETRANS
      ENVIRONMENT TRANS_USE_MPI=0 )
    if (TEST ectrans_test_transi_lam_gpu)
      set_tests_properties(ectrans_test_transi_lam_gpu PROPERTIES LABELS gpu)
    endif()
  endif()

endif()
