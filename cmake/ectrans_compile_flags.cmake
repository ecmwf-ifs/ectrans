# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


# Base flags for all supported build types
foreach( lang Fortran C CXX )
  ectrans_add_flags( "-g -O0"            BUILD DEBUG            LANG ${lang} )
  ectrans_add_flags( "-g -O2 -DNDEBUG"   BUILD RELWITHDEBINFO   LANG ${lang} )
  ectrans_add_flags( "-O3 -DNDEBUG"      BUILD RELEASE          LANG ${lang} )
  ectrans_add_flags( "-g -O2 -DNDEBUG"   BUILD BIT              LANG ${lang} )
endforeach()


# Flag to tell compiler that Fortran side has no program
# Needed if linking a C executable against some Fortran objects with some compilers
# Not needed for most
set( NO_FORTRAN_MAIN_FLAG "" )


if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  # gfortran 10 has become stricter with argument matching
  if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
    ectrans_add_fortran_flags( "-fallow-argument-mismatch" )
  endif()
endif()


if( CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
  ectrans_add_fortran_flags( "-Mlarge_arrays" )

  # should really be part of configuration, or ecbuild default?
  ectrans_add_fortran_flags( "-traceback"  BUILD DEBUG )
  ectrans_add_fortran_flags( "-fast"       BUILD RELEASE )
  ectrans_add_fortran_flags( "-gopt -fast" BUILD RELWITHDEBINFO )

  set( NO_FORTRAN_MAIN_FLAG "-Mnomain")
endif()

if( CMAKE_CXX_COMPILER_IMPORT_STD_LIBRARIES )
  # warning 177: variable x was declared but never referenced
  ectrans_add_cxx_flags( "--diag_suppress=177" )
endif()

if( CMAKE_Fortran_COMPILER_ID STREQUAL "Cray" )
  # A module named ... has already been directly or indirectly use associated into this scope
  ectrans_add_fortran_flags( "-hnomessage=878" )
  # Module ... has no public objects declared in the module, therefore nothing can be use associated
  # from the module.
  ectrans_add_fortran_flags( "-hnomessage=867" )
  # An OpenMP parallel construct in a target region is limited to a single thread.
  ectrans_add_fortran_flags( "-M7256" )
endif()


if( CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
  ectrans_add_fortran_flags( "-march=core-avx2 -no-fma" BUILD BIT )
  ectrans_add_fortran_flags( "-fp-model precise -fp-speculation=safe" )
  ectrans_add_fortran_flags( "-heap-arrays 32" )
  set( NO_FORTRAN_MAIN_FLAG "-nofor-main" )
  if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
    ectrans_add_fortran_flags( "-fast-transcendentals" )
  endif()
endif()


if( CMAKE_Fortran_COMPILER_ID STREQUAL "XL" )
  ectrans_add_fortran_flags( "-qextname -qnobindcextname" )
endif()
