#define cufftSafeCall(err) __cufftSafeCall(err, __FILE__, __LINE__)
#include "hipfft.h"
#include "stdio.h"
    static const char *_cudaGetErrorEnum(hipfftResult error)
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
    }

    return "<unknown>";
    }

    inline void __cufftSafeCall(hipfftResult err, const char *file, const int line)
    {
    if( HIPFFT_SUCCESS != err) {
    fprintf(stderr, "CUFFT error at 1\n");
    fprintf(stderr, "CUFFT error in file '%s'\n",__FILE__);
    fprintf(stderr, "CUFFT error at 2\n");
    /*fprintf(stderr, "CUFFT error line '%s'\n",__LINE__);*/
    fprintf(stderr, "CUFFT error at 3\n");
    /*fprintf(stderr, "CUFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n",__FILE__, __LINE__,err, \
    _cudaGetErrorEnum(err)); \*/
    fprintf(stderr, "CUFFT error %d: %s\nterminating!\n",err,_cudaGetErrorEnum(err)); \
    hipDeviceReset(); return; \
    }
    }

extern "C"
void
destroy_plan_fftc_(hipfftHandle *PLANp)
{
hipfftHandle plan = *PLANp;

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

cufftSafeCall(hipfftDestroy(plan));

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}


}

