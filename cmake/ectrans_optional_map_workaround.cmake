# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# CCE's offload runtime segfaults in is_contiguous_dv when a TARGET transfer list
# holds the null descriptor of an absent OPTIONAL assumed-shape dummy, which
# TRGTOL and TRLTOG hit on every call omitting a PGP* array: 138 of 366 tests.
# ECTRANS_OPTIONAL_MAP_WORKAROUND maps stand-in pointers there instead.
#
# Seen on 21.0.0; 19 and 20 want a ROCm 6 that was not available to test with.
# Set ECTRANS_HAVE_OPTIONAL_MAP_ISSUE explicitly to override.
if( NOT DEFINED ECTRANS_HAVE_OPTIONAL_MAP_ISSUE )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray"
      AND NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 20 )
    set( ECTRANS_HAVE_OPTIONAL_MAP_ISSUE True )
  else()
    set( ECTRANS_HAVE_OPTIONAL_MAP_ISSUE False )
  endif()
endif()
