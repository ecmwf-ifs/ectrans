# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( HAVE_CPU )
  # It doesn't really matter which precision we use for this test, as long as we use one that's
  # actually available
  if( HAVE_SINGLE_PRECISION )
    set( etrans_lib etrans_sp )
    set( trans_lib trans_sp )
    set( parkind_lib parkind_sp )
  else()
    set( etrans_lib etrans_dp )
    set( trans_lib trans_dp )
    set( parkind_lib parkind_dp )
  endif()
  ecbuild_add_test(
    TARGET test_etrans_trans_inq
    SOURCES test_etrans_trans_inq.F90
    LIBS ${trans_lib} ${etrans_lib} ${parkind_lib}
    LINKER_LANGUAGE Fortran
  )
endif()
