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
#include <cstdio>

// Library name
#define hipblas cublas
#define HIPBLAS CUBLAS

// CPP definitions
#define HIPBLAS_OP_T CUBLAS_OP_T
#define HIPBLAS_OP_N CUBLAS_OP_N
#define HIPBLAS_STATUS_SUCCESS CUBLAS_STATUS_SUCCESS

// Data types
#define hipError_t cudaError_t
#define hipStream_t cudaStream_t
#define hipblasHandle_t cublasHandle_t
#define hipblasStatus_t cublasStatus_t
#define hipblasOperation_t cublasOperation_t
#define hipGraph_t cudaGraph_t
#define hipGraphNode_t cudaGraphNode_t
#define hipGraphExec_t cudaGraphExec_t

// Constants
#define hipMemcpyHostToDevice cudaMemcpyHostToDevice
#define hipMemcpyDeviceToHost cudaMemcpyDeviceToHost

// Library calls
#define hipblasCreate cublasCreate
#define hipblasDestroy cublasDestroy
#define hipblasDgemm cublasDgemm
#define hipblasSgemm cublasSgemm
#define hipblasDgemmBatched cublasDgemmBatched
#define hipblasSgemmBatched cublasSgemmBatched
#define hipblasDgemmStridedBatched cublasDgemmStridedBatched
#define hipblasSgemmStridedBatched cublasSgemmStridedBatched
#define hipblasSetStream cublasSetStream

#define hipGraphExecDestroy cudaGraphExecDestroy
#define hipGraphCreate cudaGraphCreate
#define hipGraphDestroy cudaGraphDestroy
#define hipGraphLaunch cudaGraphLaunch
#define hipGraphInstantiate cudaGraphInstantiate
#define hipGraphAddChildGraphNode cudaGraphAddChildGraphNode
#define hipStreamCreate cudaStreamCreate
#define hipStreamDestroy cudaStreamDestroy
#define hipStreamCaptureModeGlobal cudaStreamCaptureModeGlobal
#define hipStreamBeginCapture cudaStreamBeginCapture
#define hipStreamEndCapture cudaStreamEndCapture

// Runtime calls
#define hipHostMalloc(PTR, SIZE, FLAGS) cudaMallocHost(PTR, SIZE)
#define hipMalloc cudaMalloc
#define hipFree cudaFree
#define hipMemcpy cudaMemcpy
#define hipDeviceSynchronize cudaDeviceSynchronize

inline static const char * _blasGetErrorEnum(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}

#define HIC_CHECK(e)                                                          \
  {                                                                            \
    cudaError_t err = (e);                                                     \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error: %s, line %d, %s: %s\n", __FILE__, __LINE__, \
              #e, cudaGetErrorString(err));                                    \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

#define HICBLAS_CHECK(e)                                                         \
{                                                                            \
  cublasStatus_t err = (e);                                                     \
  if (err != CUBLAS_STATUS_SUCCESS) {                                                  \
    fprintf(stderr, "CUDA error: %s, line %d, %s: %s\n", __FILE__, __LINE__, \
            #e, _blasGetErrorEnum(err));                                    \
    exit(EXIT_FAILURE);                                                      \
  }                                                                          \
}

#endif
