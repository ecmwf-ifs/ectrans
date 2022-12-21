#define hipfftSafeCall(err) __hipfftSafeCall(err, __FILE__, __LINE__)
#include "hip/hip_runtime.h"
#include "hipfft.h"
#include "stdio.h"
#include "execute_plan_ffth.hip.h"

#ifdef TRANS_SINGLE
typedef hipfftComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftReal HIP_DATA_TYPE_REAL;
#else
typedef hipfftDoubleComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftDoubleReal HIP_DATA_TYPE_REAL;
#endif


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

__global__ void debug(int varId, int N, HIP_DATA_TYPE_COMPLEX *x) {
    //printf("Hello from GPU\n");
    for (int i = 0; i < N; i++)
    {
        HIP_DATA_TYPE_COMPLEX a = x[i];
        double b = (double)a.x;
        double c = (double)a.y;
        if (varId == 0) printf("GPU: input[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
        if (varId == 1) printf("GPU: output[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
    }}

__global__ void debugFloat(int varId, int N, HIP_DATA_TYPE_REAL *x) {
    //printf("Hello from GPU\n");
    for (int i = 0; i < N; i++)
    {
        double a = (double)x[i];
        if (varId == 0) printf("GPU: input[%d]=%2.4f\n",i+1,a);
        if (varId == 1) printf("GPU: output[%d]=%2.4f\n",i+1,a);
    }}

/*extern "C" {

void execute_plan_ffth_c_(int ISIGNp, int N, DATA_TYPE *data_in_host, DATA_TYPE *data_out_host, long *iplan)
*/void hipfunction(int ISIGNp, int N, DATA_TYPE *data_in_host, DATA_TYPE *data_out_host, long *iplan)
{
HIP_DATA_TYPE_COMPLEX *data_in = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_in_host);
HIP_DATA_TYPE_COMPLEX *data_out = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_out_host);
hipfftHandle* PLANp = reinterpret_cast<hipfftHandle*>(iplan);
//fprintf(stderr, "execute_plan_ffth_c_: plan-address = %p\n",PLANp);
//abort();
hipfftHandle plan = *PLANp;
int ISIGN = ISIGNp;

// Check variables on the GPU:
/*int device_count = 0;
hipGetDeviceCount(&device_count);
for (int i = 0; i < device_count; ++i) {
    hipSetDevice(i);
    hipLaunchKernelGGL(debug, dim3(1), dim3(1), 0, 0, 0, N, data_in);
    hipDeviceSynchronize();
}*/

/*if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Hip error: Failed to synchronize\n");
	return;	
}*/

if( ISIGN== -1 ){
  hipfftSafeCall(hipfftExecR2C(plan, (HIP_DATA_TYPE_REAL*)data_in, data_out));
}
else if( ISIGN== 1){
  hipfftSafeCall(hipfftExecC2R(plan, data_in, (HIP_DATA_TYPE_REAL*)data_out));
}
else {
  abort();
}

hipDeviceSynchronize();

/*for (int i = 0; i < device_count; ++i) {
    hipSetDevice(i);
    hipLaunchKernelGGL(debugFloat, dim3(1), dim3(1), 0, 0, 1, N, (HIP_DATA_TYPE_REAL*)data_out);
    hipDeviceSynchronize();
}*/

//if (hipDeviceSynchronize() != hipSuccess){
//	fprintf(stderr, "Hip error: Failed to synchronize\n");
//	return;	
//}


}
//}
