// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

// Include hip runtime and rocfft headers and provide error enum translation

#ifndef __HICFFT_ROCFFT_H__
#define __HICFFT_ROCFFT_H__

#include <iostream>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#pragma-messages"
#endif
#include <hip/hip_runtime.h>
#include <rocfft/rocfft.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif

inline static const char * _fftGetErrorEnum(rocfft_status error)
{
    switch (error)
    {
        case rocfft_status_success:
        return "ROCFFT_SUCCESS";

        case rocfft_status_failure:
        return "ROCFFT_FAILURE";

        case rocfft_status_invalid_arg_value:
        return "ROCFFT_INVALID_ARG_VALUE";

        case rocfft_status_invalid_dimensions:
        return "ROCFFT_INVALID_DIMENSIONS";

        case rocfft_status_invalid_array_type:
        return "ROCFFT_INVALID_ARRAY_TYPE";

        case rocfft_status_invalid_strides:
        return "ROCFFT_INVALID_STRIDES";

        case rocfft_status_invalid_distance:
        return "ROCFFT_INVALID_DISTANCE";

        case rocfft_status_invalid_offset:
        return "ROCFFT_INVALID_OFFSET";

        case rocfft_status_invalid_work_buffer:
        return "ROCFFT_INVALID_WORK_BUFFER";
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

inline void __fftSafeCall(rocfft_status err, const char *file, const int line)
{
    if( rocfft_status_success != err) {
        fprintf(stderr, "ROCFFT error in file '%s' at line %d\n", file, line);
        fprintf(stderr, "ROCFFT error %d: %s\nterminating!\n", (int)err, _fftGetErrorEnum(err));
        std::ignore = hipDeviceReset();
        exit(EXIT_FAILURE);
    }
}

#endif

