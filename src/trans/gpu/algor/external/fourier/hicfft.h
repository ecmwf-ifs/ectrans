// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

//      HIC--->FFT
//      hip
//        cuda
//
// Common header to provide abstraction layer to utilize hipfft and cufft from
// common wrapper calls. Runtime and library specific implementations are pulled
// in from bespoke header files.
//

#ifndef __HICFFT_H__
#define __HICFFT_H__

#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include "abor1.h"

#ifdef HIPGPU
#include "hicfft_hip.h"
#elif defined(CUDAGPU)
#include "hicfft_cuda.h"
#endif

inline void _printError(const char * component, const char * file, const int line, int err, const char * err_str) {
    fprintf(stderr, "%s error at 1\n", component);
    fprintf(stderr, "%s error in file '%s'\n", component, file);
    fprintf(stderr, "%s error at 2\n", component);
    fprintf(stderr, "%s error line '%d'\n", component, line);
    fprintf(stderr, "%s error at 3\n", component);
    fprintf(stderr, "%s error %d: %s\nterminating!\n", component, err, err_str);
    return;
}

inline void __fftSafeCall(hipfftResult err, const char *file, const int line)
{
    if( hipSuccess != (int) err) {
        _printError("GPU runtime", file, line, err, _fftGetErrorEnum(err));
        hipDeviceReset();
        ABOR1("Error in FFT\n");
    }
}

#endif
