# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ectrans_find_lapack )
  # This macro sets the LAPACK_LIBRARIES variable
  # IF MKL is preferred, unless ENABLE_MKL=OFF

  if( HAVE_MKL )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
  else()
    ecbuild_find_package( NAME LAPACK REQUIRED )
  endif()
  ecbuild_debug_var( LAPACK_LIBRARIES )

endmacro()
