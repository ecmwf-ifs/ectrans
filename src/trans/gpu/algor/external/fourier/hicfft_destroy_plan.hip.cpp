#include "hicfft.h"

#define fftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)

extern "C"
void
hicfft_destroy_plan_(hipfftHandle *PLANp)
{
    hipfftHandle plan = *PLANp;

    if (hipDeviceSynchronize() != hipSuccess){
        fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
        return;
    }

    fftSafeCall(hipfftDestroy(plan));

    if (hipDeviceSynchronize() != hipSuccess){
        fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
        return;
    }

}
