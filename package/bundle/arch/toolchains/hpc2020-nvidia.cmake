# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# Compiler
####################################################################

set(CMAKE_C_COMPILER nvc)
set(CMAKE_CXX_COMPILER nvc++)
set(CMAKE_Fortran_COMPILER nvfortran)

####################################################################
# Default features
####################################################################

set( ENABLE_GPU ON )
set( ENABLE_ACC ON )
set( ENABLE_GPU_STATIC ON )
set( ENABLE_FIELD_API ON )
