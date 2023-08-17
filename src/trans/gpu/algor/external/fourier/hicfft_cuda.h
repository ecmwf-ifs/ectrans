// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

// Include cufft header and provide CPP macros to rewrite HIP and hipfft names
// to CUDA and cufft names

#ifndef __HICFFT_CUDA_H__
#define __HICFFT_CUDA_H__

#include "cufft.h"

// Library Names
#define hipfft cufft
#define HIPFFT CUFFT

// CPP macros
#define HIPFFT_SUCCESS CUFFT_SUCCESS
#define HIPFFT_R2C CUFFT_R2C
#define HIPFFT_C2R CUFFT_C2R
#define HIPFFT_D2Z CUFFT_D2Z
#define HIPFFT_Z2D CUFFT_Z2D

// Constants and types
#define hipError_t cudaError_t
#define hipSuccess cudaSuccess
#define hipfftHandle cufftHandle
#define hipfftType cufftType
#define hipfftResult cufftResult
#define hipfftComplex cufftComplex
#define hipfftReal cufftReal
#define hipfftDoubleComplex cufftDoubleComplex
#define hipfftDoubleReal cufftDoubleReal

#define hipfftCreate cufftCreate
#define hipfftDestroy cufftDestroy
#define hipfftPlanMany cufftPlanMany

#define hipfftExecR2C cufftExecR2C
#define hipfftExecC2R cufftExecC2R
#define hipfftExecD2Z cufftExecD2Z
#define hipfftExecZ2D cufftExecZ2D

// Runtime calls
#define hipDeviceSynchronize cudaDeviceSynchronize
#define hipDeviceReset cudaDeviceReset
#define _hipGetErrorEnum _cudaGetErrorEnum

inline static const char * _fftGetErrorEnum(cufftResult error)
{
    switch (error)
    {
        case CUFFT_SUCCESS:
        return "CUFFT_SUCCESS";

        case CUFFT_INVALID_PLAN:
        return "CUFFT_INVALID_PLAN";

        case CUFFT_ALLOC_FAILED:
        return "CUFFT_ALLOC_FAILED";

        case CUFFT_INVALID_TYPE:
        return "CUFFT_INVALID_TYPE";

        case CUFFT_INVALID_VALUE:
        return "CUFFT_INVALID_VALUE";

        case CUFFT_INTERNAL_ERROR:
        return "CUFFT_INTERNAL_ERROR";

        case CUFFT_EXEC_FAILED:
        return "CUFFT_EXEC_FAILED";

        case CUFFT_SETUP_FAILED:
        return "CUFFT_SETUP_FAILED";

        case CUFFT_INVALID_SIZE:
        return "CUFFT_INVALID_SIZE";

        case CUFFT_UNALIGNED_DATA:
        return "CUFFT_UNALIGNED_DATA";

        case CUFFT_INCOMPLETE_PARAMETER_LIST:
        return "CUFFT_INCOMPLETE_PARAMETER_LIST";

        case CUFFT_INVALID_DEVICE:
        return "CUFFT_INVALID_DEVICE";

        case CUFFT_PARSE_ERROR:
        return "CUFFT_PARSE_ERROR";

        case CUFFT_NO_WORKSPACE:
        return "CUFFT_NO_WORKSPACE";

        case CUFFT_NOT_IMPLEMENTED:
        return "CUFFT_NOT_IMPLEMENTED";

        case CUFFT_NOT_SUPPORTED:
        return "CUFFT_NOT_SUPPORTED";
    }

    return "<unknown>";
}

#endif
