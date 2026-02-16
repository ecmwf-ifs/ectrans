# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# CSC LUMI-G cluster, with AMD MI250X GPUs

set( OpenMP_C_LIB_NAMES craymp )
set( OpenMP_CXX_LIB_NAMES craymp )
set( OpenMP_Fortran_LIB_NAMES craymp )
set( OpenMP_craymp_LIBRARY craymp )
set( ENABLE_OMP ON )
