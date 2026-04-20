# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Check whether we have feenableexcept and fedisableexcept for controlling individual floating-point exceptions
cmake_push_check_state(RESET)
include(CheckSymbolExists)
  set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
  if(UNIX)
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif()
  check_symbol_exists(feenableexcept "fenv.h" HAVE_FEENABLEEXCEPT)
  check_symbol_exists(fedisableexcept "fenv.h" HAVE_FEDISABLEEXCEPT)
  if( HAVE_FEENABLEEXCEPT )
      set( HAVE_FEENABLEEXCEPT 1 )
  else()
      set( HAVE_FEENABLEEXCEPT 0 )
  endif()
  if( HAVE_FEDISABLEEXCEPT )
      set( HAVE_FEDISABLEEXCEPT 1 )
  else()
      set( HAVE_FEDISABLEEXCEPT 0 )
  endif()
cmake_pop_check_state()


list( APPEND util_src
    ectrans_command_line_parser_mod.F90
    ectrans_decomposition_mod.F90
    ectrans_device_mod.F90
    ectrans_error_mod.F90
    ectrans_grids_mod.F90
    ectrans_gstats_mod.F90
    ectrans_log_mod.F90
    ectrans_memory_mod.F90
    ectrans_memory.c
    ectrans_mpi_mod.F90
    ectrans_program_mod.F90
    ectrans_timer_mod.F90
    ectrans_fpe_trapping.cc
    ectrans_fpe_trapping.F90
)

list( APPEND util_PRIVATE_DEFINITIONS_GPU
  $<${HAVE_CUDA}:CUDA>
  $<${HAVE_HIP}:HIP>
)
list( APPEND util_PRIVATE_LIBRARIES_GPU
  $<${HAVE_CUDA}:CUDA::cudart>
  $<${HAVE_HIP}:hip::host>
)

ecbuild_add_library(
    TARGET  ectrans_util_cpu
    TYPE STATIC
    SOURCES ${util_src}
    PUBLIC_LIBS fiat
    PRIVATE_DEFINITIONS HAVE_FEENABLEEXCEPT=${HAVE_FEENABLEEXCEPT}
                        HAVE_FEDISABLEEXCEPT=${HAVE_FEDISABLEEXCEPT}
    LINKER_LANGUAGE Fortran
    NOINSTALL
)
ecbuild_target_fortran_module_directory(TARGET ectrans_util_cpu MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/ectrans_util_cpu )

if (HAVE_GPU)
ecbuild_add_library(
    TARGET  ectrans_util_gpu
    TYPE STATIC
    SOURCES ${util_src}
    PUBLIC_LIBS   fiat
                  $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
                  $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
                  $<${HAVE_CUDA}:CUDA::cudart>
                  $<${HAVE_HIP}:hip::host>
    PRIVATE_DEFINITIONS $<${HAVE_CUDA}:CUDA>
                        $<${HAVE_HIP}:HIP>
                        HAVE_FEENABLEEXCEPT=${HAVE_FEENABLEEXCEPT}
                        HAVE_FEDISABLEEXCEPT=${HAVE_FEDISABLEEXCEPT}
    NOINSTALL
)
ecbuild_target_fortran_module_directory(TARGET ectrans_util_gpu MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module/ectrans_util_gpu )

add_library( ectrans_util ALIAS ectrans_util_gpu )
else()
add_library( ectrans_util ALIAS ectrans_util_cpu )
endif()

