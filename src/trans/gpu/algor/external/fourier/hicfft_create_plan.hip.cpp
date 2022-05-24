#include "hicfft.h"

#define fftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)

// static int allocatedWorkspace=0;
// static void* planWorkspace;
// static int planWorkspaceSize=100*1024*1024; //100MB

extern "C"
void
hicfft_create_plan_(hipfftHandle * *plan, int *ISIGNp, int *Np, int *LOTp, int *stridep)
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

    // //create a single re-usable workspace
    // if(!allocatedWorkspace){
    //   allocatedWorkspace=1;
    //   //allocate plan workspace
    //   hipMalloc(&planWorkspace,planWorkspaceSize);
    // }
    //
    // //disable auto allocation so we can re-use a single workspace (created above)
    //  hipfftSetAutoAllocation(plan, false);

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

    if( ISIGN== -1 ){
      fftSafeCall(hipfftPlanMany(*plan, 1, &N,
                    embed, stride, dist,
                    embed, stride, dist,
                    fft_dir, LOT));
      //fftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_D2Z, LOT));
    }
    else if( ISIGN== 1){
      fftSafeCall(hipfftPlanMany(*plan, 1, &N,
                    embed, stride, dist,
                    embed, stride, dist,
                    fft_inv, LOT));
      //fftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_Z2D, LOT));
    }
    else {
      abort();
    }

    // // use our reusaable work area for the plan
    // hipfftSetWorkArea(plan,planWorkspace);

    /*
    if( ISIGN== -1 ){
      fftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_D2Z, LOT));
    }
    else if( ISIGN== 1){
      fftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_Z2D, LOT));
    }
    else {
      abort();
    }
    */

    if (hipDeviceSynchronize() != hipSuccess){
      fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
      return;
    }

    // // get size used by this plan
    // size_t workSize;
    // hipfftGetSize(plan,&workSize);
    //
    // // exit if we don't have enough space for the work area in the re-usable workspace
    // if(workSize > planWorkspaceSize){
    //   printf("create_plan_ffth: plan workspace size not large enough - exiting\n");
    // exit(1);
    // }


    return;

}
