# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI ON )
set( ENABLE_ACC OFF )

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffpe-trap=invalid,zero,overflow")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fstack-arrays")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fconvert=big-endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fbacktrace")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fno-second-underscore")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffree-form")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffast-math")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fno-unsafe-math-optimizations")
