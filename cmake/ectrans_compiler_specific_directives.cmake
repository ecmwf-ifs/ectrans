# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Per-compiler spellings for the OpenMP clauses in the GPU backend. The defaults are the
# amdflang spellings the directives were written against; nvfortran needs all four replaced,
# Cray Fortran only the DEFAULT argument:
#
#   ECTRANS_MAP_PRESENT_ALLOC   Map-type modifiers for MAP(...) of device-resident storage.
#                               nvfortran has no OpenMP 5.1 'present' modifier (absent in 26.3
#                               and 26.5); it is only a runtime check, so dropping it leaves
#                               the mapping semantics unchanged.
#   ECTRANS_OMP_DEFAULT         Argument of DEFAULT(...). nvfortran's DEFAULT(NONE) check
#                               cannot resolve the ASSOCIATE selectors these constructs use
#                               (R_NTMAX=>R%NTMAX and similar), whatever SHARED lists;
#                               DEFAULT(SHARED) is what the implicit rules give anyway.
#                               Cray Fortran 21 rejects the same ASSOCIATE selectors, and in
#                               addition wants an explicit scope for the HAS_DEVICE_ADDR list
#                               items of TRLTOM and TRGTOL. Naming those in SHARED() is not a
#                               way out either: the UPDSPB construct already carries a comment
#                               that a longer SHARED list there hits ftn-7991.
#   ECTRANS_DEVICE_ADDR_CLAUSE  Clause naming storage that is already on the device.
#                               nvfortran implements neither HAS_DEVICE_ADDR nor a usable
#                               substitute: IS_DEVICE_PTR parses but carries no array
#                               descriptor, so assumed-shape list items fault. SHARED is
#                               correct but copies the array in and out per launch wherever
#                               the runtime cannot resolve the storage; LPGP_ON_GPU removes
#                               that traffic for the gridpoint arrays.
#   ECTRANS_LOOP_BOUNDS_CLAUSE  Clause for read-only loop-bound scalars. On nvfortran this
#                               decides the launch geometry: FIRSTPRIVATE emits a grid-stride
#                               kernel of one block per SM however large the nest, SHARED a
#                               trip-count-sized grid. Only read-only bounds may move, as a
#                               scalar in SHARED is still implicitly firstprivate on the
#                               target construct and so needs no MAP.
if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  set( ECTRANS_MAP_PRESENT_ALLOC  "ALLOC" )
  set( ECTRANS_OMP_DEFAULT        "SHARED" )
  set( ECTRANS_DEVICE_ADDR_CLAUSE "SHARED" )
  set( ECTRANS_LOOP_BOUNDS_CLAUSE "SHARED" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( ECTRANS_MAP_PRESENT_ALLOC  "PRESENT,ALLOC" )
  set( ECTRANS_OMP_DEFAULT        "SHARED" )
  set( ECTRANS_DEVICE_ADDR_CLAUSE "HAS_DEVICE_ADDR" )
  set( ECTRANS_LOOP_BOUNDS_CLAUSE "FIRSTPRIVATE" )
else()
  set( ECTRANS_MAP_PRESENT_ALLOC  "PRESENT,ALLOC" )
  set( ECTRANS_OMP_DEFAULT        "NONE" )
  set( ECTRANS_DEVICE_ADDR_CLAUSE "HAS_DEVICE_ADDR" )
  set( ECTRANS_LOOP_BOUNDS_CLAUSE "FIRSTPRIVATE" )
endif()

# ECTRANS_OMP_DEFAULT_CLAUSE is the whole DEFAULT clause rather than its argument, so that it
# can expand to nothing: amdflang 23.3.0 and 24.1.0-pre reject DEFAULT alongside
# HAS_DEVICE_ADDR. Being possibly empty it must share a line with other clauses rather than sit
# on a continuation of its own. Constructs without ECTRANS_DEVICE_ADDR_CLAUSE keep plain
# DEFAULT(ECTRANS_OMP_DEFAULT).
if( CMAKE_Fortran_COMPILER_ID MATCHES "Flang" )
  set( ECTRANS_OMP_DEFAULT_CLAUSE "" )
else()
  set( ECTRANS_OMP_DEFAULT_CLAUSE "DEFAULT(${ECTRANS_OMP_DEFAULT})" )
endif()

# The PGP* gridpoint arrays of TRGTOL and TRLTOG are deliberately not covered by either
# macro above. They take a plain MAP(ALLOC:...) on every compiler, which is neither of the
# two expansions here, for reasons written out at the constructs that use them.
