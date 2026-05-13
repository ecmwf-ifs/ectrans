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

if( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  ecbuild_add_fortran_flags( "-qextname -qnobindcextname" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  # gfortran 10 has become stricter with argument matching
  if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
    ecbuild_add_fortran_flags( "-fallow-argument-mismatch" )
  endif()
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  ecbuild_add_fortran_flags( "-Mlarge_arrays" )

  # should really be part of configuration, or ecbuild default?
  ecbuild_add_fortran_flags( "-traceback"      BUILD DEBUG )
  ecbuild_add_fortran_flags( "-fast"           BUILD RELEASE )
  ecbuild_add_fortran_flags( "-gopt -fast"     BUILD RELWITHDEBINFO )

  set( NO_FORTRAN_MAIN_FLAG "-Mnomain")
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  # A module named ... has already been directly or indirectly use associated into this scope
  ecbuild_add_fortran_flags( "-hnomessage=878" )
  # Module ... has no public objects declared in the module, therefore nothing can be use associated
  # from the module.
  ecbuild_add_fortran_flags( "-hnomessage=867" )
  # An OpenMP parallel construct in a target region is limited to a single thread.
  ecbuild_add_fortran_flags( "-M7256" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM" )
  ecbuild_add_fortran_flags( "-march=core-avx2 -no-fma" BUILD BIT )
  ecbuild_add_fortran_flags( "-fp-model precise -fp-speculation=safe" )
  set( NO_FORTRAN_MAIN_FLAG "-nofor-main" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  ecbuild_add_fortran_flags( "-march=core-avx2 -no-fma" BUILD BIT )
  ecbuild_add_fortran_flags( "-fast-transcendentals -fp-model precise -fp-speculation=safe" )
  set( NO_FORTRAN_MAIN_FLAG "-nofor-main" )
endif()

