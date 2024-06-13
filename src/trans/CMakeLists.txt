# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( HAVE_CPU)
  add_subdirectory( cpu )
endif()
if( HAVE_GPU )
  add_subdirectory( gpu )
endif()

## Install trans interface

file( GLOB trans_interface include/ectrans/* )
install(
  FILES        ${trans_interface}
  DESTINATION  include/ectrans
)
