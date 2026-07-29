# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Write CMake-resolved dependency values for the CI provenance collector. These
# values are normal CMake variables, so they are not necessarily in CMakeCache.

file(GENERATE
  OUTPUT "${CMAKE_BINARY_DIR}/ectrans-reference-provenance.txt"
  CONTENT
"LAPACK_LIBRARIES=${LAPACK_LIBRARIES}
LAPACK_VERSION=${LAPACK_VERSION}
FFTW_LIBRARIES=${FFTW_LIBRARIES}
FFTW_VERSION=${FFTW_VERSION}
ECTRANS_GPU_HIP_LIBRARIES=${ECTRANS_GPU_HIP_LIBRARIES}
CUDA_CUBLAS_LIBRARY=${CUDA_cublas_LIBRARY}
CUDA_CUFFT_LIBRARY=${CUDA_cufft_LIBRARY}
CMAKE_VERSION=${CMAKE_VERSION}
CMAKE_C_COMPILER_ID=${CMAKE_C_COMPILER_ID}
CMAKE_C_COMPILER_VERSION=${CMAKE_C_COMPILER_VERSION}
CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}
CMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}
CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}
CMAKE_CUDA_COMPILER_ID=${CMAKE_CUDA_COMPILER_ID}
CMAKE_CUDA_COMPILER_VERSION=${CMAKE_CUDA_COMPILER_VERSION}
CMAKE_HIP_COMPILER_ID=${CMAKE_HIP_COMPILER_ID}
CMAKE_HIP_COMPILER_VERSION=${CMAKE_HIP_COMPILER_VERSION}
ECBUILD_VERSION=${ecbuild_VERSION}
FIAT_VERSION=${fiat_VERSION}
FIELD_API_VERSION=${field_api_VERSION}
")
