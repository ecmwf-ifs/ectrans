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