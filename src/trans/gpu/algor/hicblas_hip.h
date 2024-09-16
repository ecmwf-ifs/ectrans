// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

// Include hip runtime and hipblas headers

#ifndef __HICBLAS_HIP_H__
#define __HICBLAS_HIP_H__

#include <cstdio>
#include <hip/hip_runtime.h>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#pragma-messages"
#endif
#include "hipblas/hipblas.h"
#ifdef __clang__
#pragma clang diagnostic pop
#endif

inline static const char * _blasGetErrorEnum(hipblasStatus_t error)
{
    switch (error)
    {
       case HIPBLAS_STATUS_SUCCESS:
            return "HIPBLAS_STATUS_SUCCESS";

        case HIPBLAS_STATUS_NOT_INITIALIZED:
            return "HIPBLAS_STATUS_NOT_INITIALIZED";

        case HIPBLAS_STATUS_ALLOC_FAILED:
            return "HIPBLAS_STATUS_ALLOC_FAILED";

        case HIPBLAS_STATUS_INVALID_VALUE:
            return "HIPBLAS_STATUS_INVALID_VALUE";

        case HIPBLAS_STATUS_ARCH_MISMATCH:
            return "HIPBLAS_STATUS_ARCH_MISMATCH";

        case HIPBLAS_STATUS_MAPPING_ERROR:
            return "HIPBLAS_STATUS_MAPPING_ERROR";

        case HIPBLAS_STATUS_EXECUTION_FAILED:
            return "HIPBLAS_STATUS_EXECUTION_FAILED";

        case HIPBLAS_STATUS_INTERNAL_ERROR:
            return "HIPBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}

#define HICBLAS_CHECK(e)                                                      \
  {                                                                           \
    hipblasStatus_t err = (e);                                                \
    if (err != HIPBLAS_STATUS_SUCCESS) {                                      \
      fprintf(stderr, "HIP error: %s, line %d, %s: %s\n", __FILE__, __LINE__, \
              #e, _blasGetErrorEnum(err));                                    \
      exit(EXIT_FAILURE);                                                     \
    }                                                                         \
  }

#define HIC_CHECK(e)                                                          \
  {                                                                           \
    hipError_t err = (e);                                                     \
    if (err != hipSuccess) {                                                  \
      fprintf(stderr, "HIP error: %s, line %d, %s: %s\n", __FILE__, __LINE__, \
              #e, hipGetErrorString(err));                                    \
      exit(EXIT_FAILURE);                                                     \
    }                                                                         \
  }

#endif
