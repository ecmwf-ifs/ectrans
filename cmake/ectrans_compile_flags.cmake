# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Flag to tell compiler that Fortran side has no program
# Needed if linking a C executable against some Fortran objects with some compilers
# Not needed for most
set( NO_FORTRAN_MAIN_FLAG "" )

if( CMAKE_Fortran_COMPILER_ID STREQUAL "XL" )
  ecbuild_add_fortran_flags( "-qextname -qnobindcextname" )
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  # gfortran 10 has become stricter with argument matching
  if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
    ecbuild_add_fortran_flags( "-fallow-argument-mismatch" )
  endif()
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
  ecbuild_add_fortran_flags( "-Mlarge_arrays" )

  # should really be part of configuration, or ecbuild default?
  ecbuild_add_fortran_flags( "-traceback"      BUILD DEBUG )
  ecbuild_add_fortran_flags( "-fast"           BUILD RELEASE )
  ecbuild_add_fortran_flags( "-gopt -fast"     BUILD RELWITHDEBINFO )

  set( NO_FORTRAN_MAIN_FLAG "-Mnomain")
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Cray" )
  # A module named ... has already been directly or indirectly use associated into this scope
  ecbuild_add_fortran_flags( "-hnomessage=878" )
  # Module ... has no public objects declared in the module, therefore nothing can be use associated
  # from the module.
  ecbuild_add_fortran_flags( "-hnomessage=867" )
  # An OpenMP parallel construct in a target region is limited to a single thread.
  ecbuild_add_fortran_flags( "-M7256" )
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
  ecbuild_add_fortran_flags( "-march=core-avx2 -no-fma" BUILD BIT )
  ecbuild_add_fortran_flags( "-fp-model precise -fp-speculation=safe" )
  ecbuild_add_fortran_flags( "-heap-arrays 32" )
  set( NO_FORTRAN_MAIN_FLAG "-nofor-main" )
  if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
    ecbuild_add_fortran_flags( "-fast-transcendentals" )
  endif()
endif()

ecbuild_add_fortran_flags( "-O3 -DNDEBUG" NAME base_release BUILD RELEASE )
ecbuild_add_c_flags( "-O3 -DNDEBUG" NAME base_release BUILD RELEASE )
ecbuild_add_cxx_flags( "-O3 -DNDEBUG" NAME base_release BUILD RELEASE )

ecbuild_add_fortran_flags( "-O2 -DNDEBUG" NAME base_bit BUILD BIT )
ecbuild_add_c_flags( "-O2 -DNDEBUG" NAME base_bit BUILD BIT )
ecbuild_add_cxx_flags( "-O2 -DNDEBUG" NAME base_bit BUILD BIT )

ecbuild_add_fortran_flags( "-g -O0"   NAME base_debug BUILD DEBUG )
ecbuild_add_c_flags( "-g -O0"   NAME base_debug BUILD DEBUG )
ecbuild_add_cxx_flags( "-g -O0"   NAME base_debug BUILD DEBUG )

ecbuild_add_fortran_flags( "-g -O2 -DNDEBUG" NAME base_relwithdebinfo BUILD RELWITHDEBINFO )
ecbuild_add_c_flags( "-g -O2 -DNDEBUG" NAME base_relwithdebinfo BUILD RELWITHDEBINFO )
ecbuild_add_cxx_flags( "-g -O2 -DNDEBUG" NAME base_relwithdebinfo BUILD RELWITHDEBINFO )
