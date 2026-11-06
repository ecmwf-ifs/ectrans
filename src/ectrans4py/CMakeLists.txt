# (C) Copyright 2025- ECMWF.
# (C) Copyright 2025- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# (using CMAKE_CURRENT_SOURCE_DIR is necessary because sources are in a different directory than the
# target library (trans_${prec}))
ecbuild_list_add_pattern(
  LIST ectrans4py_src
  GLOB ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
  QUIET
)

if( HAVE_DOUBLE_PRECISION )
  set( prec "dp" )
  set( prec_cpp "DOUBLE_PRECISION")
else()
  set( prec "sp" )
  set( prec_cpp "SINGLE_PRECISION")
endif()

# Add ectrans4py library
ecbuild_add_library(
      TARGET           ectrans4py_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${ectrans4py_src}
      PUBLIC_LIBS      fiat trans_${prec} parkind_${prec} etrans_${prec}
)
target_compile_definitions( ectrans4py_${prec} PRIVATE ${prec_cpp} )

# OpenMP threads the blocked repacking of spectral arrays (spec_layout4py.F90, with the
# environment variable ECTRANS4PY_OMP set). The OpenMP compile flag goes on that source only, the
# only one with OpenMP directives, and only the OpenMP libraries on the target: with the flag on
# every source (which linking OpenMP::OpenMP_Fortran would do), the code generated for the other
# sources changes: with the Intel compiler, that made one of the transforms tens of times slower.
if( HAVE_OMP )
  separate_arguments( ectrans4py_omp_flags NATIVE_COMMAND "${OpenMP_Fortran_FLAGS}" )
  # (ECTRANS4PY_WITH_OPENMP: the build asks for OpenMP, which the source checks it is compiled with)
  set_source_files_properties( ${CMAKE_CURRENT_SOURCE_DIR}/spec_layout4py.F90
    PROPERTIES COMPILE_OPTIONS "${ectrans4py_omp_flags}"
               COMPILE_DEFINITIONS ECTRANS4PY_WITH_OPENMP )
  if( OpenMP_Fortran_LIBRARIES )
    target_link_libraries( ectrans4py_${prec} PRIVATE ${OpenMP_Fortran_LIBRARIES} )
  else()
    target_link_options( ectrans4py_${prec} PRIVATE ${ectrans4py_omp_flags} )
  endif()
endif()