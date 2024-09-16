// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

// Include hip runtime and hipfft headers and provide error enum translation

#ifndef __HICFFT_HIP_H__
#define __HICFFT_HIP_H__

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#pragma-messages"
#endif
#include <hip/hip_runtime.h>
#include "hipfft/hipfft.h"
#ifdef __clang__
#pragma clang diagnostic pop
#endif

inline static const char * _fftGetErrorEnum(hipfftResult error)
{
    switch (error)
    {
        case HIPFFT_SUCCESS:
        return "HIPFFT_SUCCESS";

        case HIPFFT_INVALID_PLAN:
        return "HIPFFT_INVALID_PLAN";

        case HIPFFT_ALLOC_FAILED:
        return "HIPFFT_ALLOC_FAILED";

        case HIPFFT_INVALID_TYPE:
        return "HIPFFT_INVALID_TYPE";

        case HIPFFT_INVALID_VALUE:
        return "HIPFFT_INVALID_VALUE";

        case HIPFFT_INTERNAL_ERROR:
        return "HIPFFT_INTERNAL_ERROR";

        case HIPFFT_EXEC_FAILED:
        return "HIPFFT_EXEC_FAILED";

        case HIPFFT_SETUP_FAILED:
        return "HIPFFT_SETUP_FAILED";

        case HIPFFT_INVALID_SIZE:
        return "HIPFFT_INVALID_SIZE";

        case HIPFFT_UNALIGNED_DATA:
        return "HIPFFT_UNALIGNED_DATA";

        case HIPFFT_INCOMPLETE_PARAMETER_LIST:
        return "HIPFFT_INCOMPLETE_PARAMETER_LIST";

        case HIPFFT_INVALID_DEVICE:
        return "HIPFFT_INVALID_DEVICE";

        case HIPFFT_PARSE_ERROR:
        return "HIPFFT_PARSE_ERROR";

        case HIPFFT_NO_WORKSPACE:
        return "HIPFFT_NO_WORKSPACE";

        case HIPFFT_NOT_IMPLEMENTED:
        return "HIPFFT_NOT_IMPLEMENTED";

        case HIPFFT_NOT_SUPPORTED:
        return "HIPFFT_NOT_SUPPORTED";
    }

    return "<unknown>";
}
  
#define HIC_CHECK(e)                                                         \
  {                                                                            \
    hipError_t err = (e);                                                     \
    if (err != hipSuccess) {                                                  \
      fprintf(stderr, "HIP error: %s, line %d, %s: %s\n", __FILE__, __LINE__, \
              #e, hipGetErrorString(err));                                    \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }


#endif
