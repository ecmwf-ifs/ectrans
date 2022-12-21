#define hipfftSafeCall(err) __hipfftSafeCall(err, __FILE__, __LINE__)
#include "hipfft.h"
#include "stdio.h"
    static const char *_hipGetErrorEnum(hipfftResult error)
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

    inline void __hipfftSafeCall(hipfftResult err, const char *file, const int line)
    {
    if( HIPFFT_SUCCESS != err) {
    fprintf(stderr, "HIPFFT error at 1\n");
    fprintf(stderr, "HIPFFT error in file '%s'\n",__FILE__);
    fprintf(stderr, "HIPFFT error at 2\n");
    /*fprintf(stderr, "HIPFFT error line '%s'\n",__LINE__);*/
    fprintf(stderr, "HIPFFT error at 3\n");
    /*fprintf(stderr, "HIPFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n",__FILE__, __LINE__,err, \
    _hipGetErrorEnum(err)); \*/
    fprintf(stderr, "HIPFFT error %d: %s\nterminating!\n",err,_hipGetErrorEnum(err)); \
    hipDeviceReset(); return; \
    }
    }

extern "C"
void
destroy_plan_ffth_(hipfftHandle *PLANp)
{
hipfftHandle plan = *PLANp;

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Hip error: Failed to synchronize\n");
	return;	
}

hipfftSafeCall(hipfftDestroy(plan));

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Hip error: Failed to synchronize\n");
	return;	
}


}

