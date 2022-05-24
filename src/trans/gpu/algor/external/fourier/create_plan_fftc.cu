#define cufftSafeCall(err) __cufftSafeCall(err, __FILE__, __LINE__)
#include "cufft.h"
#include "stdio.h"
    static const char *_cudaGetErrorEnum(cufftResult error)
    {
    switch (error)
    {
    case CUFFT_SUCCESS:
    return "CUFFT_SUCCESS";

    case CUFFT_INVALID_PLAN:
    return "CUFFT_INVALID_PLAN";

    case CUFFT_ALLOC_FAILED:
    return "CUFFT_ALLOC_FAILED";

    case CUFFT_INVALID_TYPE:
    return "CUFFT_INVALID_TYPE";

    case CUFFT_INVALID_VALUE:
    return "CUFFT_INVALID_VALUE";

    case CUFFT_INTERNAL_ERROR:
    return "CUFFT_INTERNAL_ERROR";

    case CUFFT_EXEC_FAILED:
    return "CUFFT_EXEC_FAILED";

    case CUFFT_SETUP_FAILED:
    return "CUFFT_SETUP_FAILED";

    case CUFFT_INVALID_SIZE:
    return "CUFFT_INVALID_SIZE";

    case CUFFT_UNALIGNED_DATA:
    return "CUFFT_UNALIGNED_DATA";
    }

    return "<unknown>";
    }

    inline void __cufftSafeCall(cufftResult err, const char *file, const int line)
    {
    if( CUFFT_SUCCESS != err) {
    fprintf(stderr, "CUFFT error at 1\n");
    fprintf(stderr, "CUFFT error in file '%s'\n",__FILE__);
    fprintf(stderr, "CUFFT error at 2\n");
    /*fprintf(stderr, "CUFFT error line '%s'\n",__LINE__);*/
    fprintf(stderr, "CUFFT error at 3\n");
    /*fprintf(stderr, "CUFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n",__FILE__, __LINE__,err, \
    _cudaGetErrorEnum(err)); \*/
    fprintf(stderr, "CUFFT error %d: %s\nterminating!\n",err,_cudaGetErrorEnum(err)); \
    cudaDeviceReset(); return; \
    }
    }


static int allocatedWorkspace=0;
static void* planWorkspace;
static int planWorkspaceSize=100*1024*1024; //100MB
 
extern "C"
void
create_plan_fftc_(cufftHandle *PLANp, int *ISIGNp, int *Np, int *LOTp, int *stridep)
{
int ISIGN = *ISIGNp;
int N = *Np;
int LOT = *LOTp;
int stride = *stridep;

cufftHandle plan;

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}


// //create a single re-usable workspace
// if(!allocatedWorkspace){
//   allocatedWorkspace=1;
//   //allocate plan workspace
//   cudaMalloc(&planWorkspace,planWorkspaceSize);
// }
//
// //disable auto allocation so we can re-use a single workspace (created above)
//  cufftSetAutoAllocation(plan, false);

int embed[1];
int dist;

#ifdef TRANS_SINGLE
cufftType cufft_1 = CUFFT_R2C;
cufftType cufft_2 = CUFFT_C2R;
#else
cufftType cufft_1 = CUFFT_D2Z;
cufftType cufft_2 = CUFFT_Z2D;
#endif

embed[0] = 1;
dist     = 1;

cufftSafeCall(cufftCreate(&plan));

//printf("CreatePlan cuFFT\n","N=",N);
//printf("%s %d \n","plan=",plan);
//printf("%s %d \n","LOT=",LOT);
//printf("%s %d \n","ISIGN=",ISIGN);
//printf("%s %d \n","Np=",*Np);

if( ISIGN== -1 ){
  cufftSafeCall(cufftPlanMany(&plan, 1, &N,
                 embed, stride, dist, 
                 embed, stride, dist, 
                 cufft_1, LOT));
  //cufftSafeCall(cufftPlan1d(&plan, N, CUFFT_D2Z, LOT));
}
else if( ISIGN== 1){
  cufftSafeCall(cufftPlanMany(&plan, 1, &N,
                 embed, stride, dist, 
                 embed, stride, dist, 
                 cufft_2, LOT));
  //cufftSafeCall(cufftPlan1d(&plan, N, CUFFT_Z2D, LOT));
}
else {
  abort();
}

// // use our reusaable work area for the plan
// cufftSetWorkArea(plan,planWorkspace); 

/*
if( ISIGN== -1 ){
  cufftSafeCall(cufftPlan1d(&plan, N, CUFFT_D2Z, LOT));
}
else if( ISIGN== 1){
  cufftSafeCall(cufftPlan1d(&plan, N, CUFFT_Z2D, LOT));
}
else {
  abort();
}
*/

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

*PLANp=plan;

// // get size used by this plan
// size_t workSize;
// cufftGetSize(plan,&workSize);
//
// // exit if we don't have enough space for the work area in the re-usable workspace
// if(workSize > planWorkspaceSize){
//   printf("create_plan_fftc: plan workspace size not large enough - exiting\n");
// exit(1);
// }


return;


}

