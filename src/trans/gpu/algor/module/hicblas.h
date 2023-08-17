// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

//      HIC--->BLAS
//      hip
//        cuda
//
// Common header to provide abstraction layer to utilize hipblas and cublas from
// common wrapper calls. Runtime and library specific implementations are pulled
// in from bespoke header files.
//

#ifndef __HICBLAS_H__
#define __HICBLAS_H__

#ifdef HIPGPU
#include "hicblas_hip.h"
#elif defined(CUDAGPU)
#include "hicblas_cuda.h"
#endif

#endif
