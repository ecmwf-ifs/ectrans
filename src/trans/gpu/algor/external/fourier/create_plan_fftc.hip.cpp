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


static int allocatedWorkspace=0;
static void* planWorkspace;
static int planWorkspaceSize=100*1024*1024; //100MB
 
extern "C"
void
create_plan_fftc_(hipfftHandle *PLANp, int *ISIGNp, int *Np, int *LOTp)
{
int ISIGN = *ISIGNp;
int N = *Np;
int LOT = *LOTp;

hipfftHandle plan;

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
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
int stride;
int dist;

#ifdef TRANS_SINGLE
hipfftType cufft_1 = HIPFFT_R2C;
hipfftType cufft_2 = HIPFFT_C2R;
#else
hipfftType cufft_1 = HIPFFT_D2Z;
hipfftType cufft_2 = HIPFFT_Z2D;
#endif

embed[0] = 1;
stride   = LOT;
dist     = 1;

cufftSafeCall(hipfftCreate(&plan));

//printf("CreatePlan cuFFT\n","N=",N);
//printf("%s %d \n","plan=",plan);
//printf("%s %d \n","LOT=",LOT);
//printf("%s %d \n","ISIGN=",ISIGN);
//printf("%s %d \n","Np=",*Np);

if( ISIGN== -1 ){
  cufftSafeCall(hipfftPlanMany(&plan, 1, &N,
                 embed, stride, dist, 
                 embed, stride, dist, 
                 cufft_1, LOT));
  //cufftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_D2Z, LOT));
}
else if( ISIGN== 1){
  cufftSafeCall(hipfftPlanMany(&plan, 1, &N,
                 embed, stride, dist, 
                 embed, stride, dist, 
                 cufft_2, LOT));
  //cufftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_Z2D, LOT));
}
else {
  abort();
}

// // use our reusaable work area for the plan
// hipfftSetWorkArea(plan,planWorkspace); 

/*
if( ISIGN== -1 ){
  cufftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_D2Z, LOT));
}
else if( ISIGN== 1){
  cufftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_Z2D, LOT));
}
else {
  abort();
}
*/

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

*PLANp=plan;

// // get size used by this plan
// size_t workSize;
// hipfftGetSize(plan,&workSize);
//
// // exit if we don't have enough space for the work area in the re-usable workspace
// if(workSize > planWorkspaceSize){
//   printf("create_plan_fftc: plan workspace size not large enough - exiting\n");
// exit(1);
// }


return;


}

