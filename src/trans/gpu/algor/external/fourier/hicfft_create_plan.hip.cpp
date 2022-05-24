#include "hicfft.h"

#define fftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)

// static int allocatedWorkspace=0;
// static void* planWorkspace;
// static int planWorkspaceSize=100*1024*1024; //100MB
void *planWorkspace;
static int currentWorkspaceSize = 0;

extern "C"
void
hicfft_create_plan_(hipfftHandle * *plan, int *ISIGNp, int *Np, int *LOTp, int *stridep, int *plan_size)
{
    int ISIGN = *ISIGNp;
    int N = *Np;
    int LOT = *LOTp;
    int stride = *stridep;

    *plan = new hipfftHandle;

    if (hipDeviceSynchronize() != hipSuccess){
      fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
      return;
    }

    int embed[1];
    int dist;

    #ifdef TRANS_SINGLE
    hipfftType fft_dir = HIPFFT_R2C;
    hipfftType fft_inv = HIPFFT_C2R;
    #else
    hipfftType fft_dir = HIPFFT_D2Z;
    hipfftType fft_inv = HIPFFT_Z2D;
    #endif

    embed[0] = 1;
    dist     = 1;

    fftSafeCall(hipfftCreate(*plan));

    // Disable auto allocation
    fftSafeCall(hipfftSetAutoAllocation(**plan, false));

    if( ISIGN== -1 ){
      fftSafeCall(hipfftPlanMany(*plan, 1, &N,
                    embed, stride, dist,
                    embed, stride/2, dist,
                    fft_dir, LOT));
    }
    else if( ISIGN== 1){
      fftSafeCall(hipfftPlanMany(*plan, 1, &N,
                    embed, stride/2, dist,
                    embed, stride, dist,
                    fft_inv, LOT));
    }
    else {
      abort();
    }

    // get size used by this plan
    size_t thisWorkplanSize;
    hipfftGetSize(**plan, &thisWorkplanSize);

    // check if this the work space is sufficiently large
    if (thisWorkplanSize > currentWorkspaceSize) {
      hipDeviceSynchronize();
      hipFree(planWorkspace);
      hipMalloc(&planWorkspace, thisWorkplanSize);
      currentWorkspaceSize = thisWorkplanSize;
    }

    if (hipDeviceSynchronize() != hipSuccess){
      fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
      return;
    }

    return;

}
