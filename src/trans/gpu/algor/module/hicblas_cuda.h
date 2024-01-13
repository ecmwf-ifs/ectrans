// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

// Include cublas header and provide CPP macros to rewrite HIP and hipblas names
// to CUDA and cublas names

#ifndef __HICBLAS_CUDA_H__
#define __HICBLAS_CUDA_H__

#include "cublas_v2.h"

// Library name
#define hipblas cublas
#define HIPBLAS CUBLAS

// CPP definitions
#define HIPBLAS_OP_T CUBLAS_OP_T
#define HIPBLAS_OP_N CUBLAS_OP_N
#define HIPBLAS_STATUS_SUCCESS CUBLAS_STATUS_SUCCESS

// Data types
#define hipError_t cudaError_t
#define hipblasHandle_t cublasHandle_t
#define hipblasStatus_t cublasStatus_t
#define hipblasOperation_t cublasOperation_t

// Constants
#define hipMemcpyHostToDevice cudaMemcpyHostToDevice
#define hipMemcpyDeviceToHost cudaMemcpyDeviceToHost

// Library calls
#define hipblasCreate cublasCreate
#define hipblasDestroy cublasDestroy
#define hipblasDgemmBatched cublasDgemmBatched
#define hipblasSgemmBatched cublasSgemmBatched
#define hipblasDgemmStridedBatched cublasDgemmStridedBatched
#define hipblasSgemmStridedBatched cublasSgemmStridedBatched

// Runtime calls
#define hipHostMalloc(PTR, SIZE, FLAGS) cudaMallocHost(PTR, SIZE)
#define hipMalloc cudaMalloc
#define hipFree cudaFree
#define hipMemcpy cudaMemcpy
#define hipDeviceSynchronize cudaDeviceSynchronize

#endif
