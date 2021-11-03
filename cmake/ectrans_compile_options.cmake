# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


if( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  ecbuild_add_fortran_flags("-qextname -qnobindcextname")
endif()

# gfortran 10 has become stricter with argument matching
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU"
    AND NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
  ecbuild_add_fortran_flags("-fallow-argument-mismatch")
endif()
