# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# See https://github.com/ecmwf-ifs/ectrans/pull/98
# There is a problem with CONTIGUOUS keyword in dist_spec_control_mod.F90
if( NOT DEFINED ECTRANS_HAVE_CONTIGUOUS_ISSUE )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel"  )
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS_EQUAL 19)
      set( ECTRANS_HAVE_CONTIGUOUS_ISSUE True )
    endif()
  elseif( CMAKE_Fortran_COMPILER_ID MATCHES "GNU"  )
    # GCC versions 9.2, 11.2, 12.2, 13.3, 14.2 are all known to have an issue with `contiguous`
    # Logic below is defensive and assumes future versions of gcc are likely to also have the issue
    if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 9 )
      set( ECTRANS_HAVE_CONTIGUOUS_ISSUE True )
    endif()
  endif()
endif()
