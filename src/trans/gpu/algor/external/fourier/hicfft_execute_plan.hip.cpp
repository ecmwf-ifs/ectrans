#include "hicfft.h"

#define fftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)

#ifdef TRANS_SINGLE
typedef float DATA_TYPE;
typedef hipfftComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftReal HIP_DATA_TYPE_REAL;
#define fftExecDir hipfftExecR2C
#define fftExecInv hipfftExecC2R
#else
typedef double DATA_TYPE;
typedef hipfftDoubleComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftDoubleReal HIP_DATA_TYPE_REAL;
#define fftExecDir hipfftExecD2Z
#define fftExecInv hipfftExecZ2D
#endif

__global__ void debug(int varId, int N, HIP_DATA_TYPE_COMPLEX *x) {
    for (int i = 0; i < N; i++)
    {
        HIP_DATA_TYPE_COMPLEX a = x[i];
        double b = (double)a.x;
        double c = (double)a.y;
        if (varId == 0) printf("GPU: input[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
        if (varId == 1) printf("GPU: output[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
    }
}

__global__ void debugFloat(int varId, int N, HIP_DATA_TYPE_REAL *x) {
    for (int i = 0; i < N; i++)
    {
        double a = (double)x[i];
        if (varId == 0) printf("GPU: input[%d]=%2.4f\n",i+1,a);
        if (varId == 1) printf("GPU: output[%d]=%2.4f\n",i+1,a);
    }
}

extern "C"
void
hicfft_execute_plan_(int ISIGNp, int N, DATA_TYPE *data_in_host, DATA_TYPE *data_out_host, long *iplan)
{
    HIP_DATA_TYPE_COMPLEX *data_in = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_in_host);
    HIP_DATA_TYPE_COMPLEX *data_out = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_out_host);
    hipfftHandle* PLANp = reinterpret_cast<hipfftHandle*>(iplan);
    hipfftHandle plan = *PLANp;
    int ISIGN = ISIGNp;

    /*if (hipDeviceSynchronize() != hipSuccess){
        fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
        return;
    }*/

    if( ISIGN== -1 ){
        fftSafeCall(fftExecDir(plan, (HIP_DATA_TYPE_REAL*)data_in, data_out));
    }
    else if( ISIGN== 1){
        fftSafeCall(fftExecInv(plan, data_in, (HIP_DATA_TYPE_REAL*)data_out));
    }
    else {
        abort();
    }

    if (hipDeviceSynchronize() != hipSuccess){
        fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
        return;
    }

}
