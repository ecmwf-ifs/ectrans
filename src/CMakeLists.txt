# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

add_subdirectory( trans )
add_subdirectory( programs )
if( HAVE_TRANSI )
  add_subdirectory(transi)
endif()
if( HAVE_ETRANS )
  add_subdirectory(etrans)
endif()
if(HAVE_ECTRANS4PY)
  add_subdirectory(ectrans4py)
endif()
