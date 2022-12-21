#define hipfftSafeCall(err) __hipfftSafeCall(err, __FILE__, __LINE__)
#include "stdio.h"
#include <hip/hip_runtime.h>
#include "hipfft.h"
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


static int allocatedWorkspace=0;
static void* planWorkspace;
static int planWorkspaceSize=100*1024*1024; //100MB
 
extern "C"
void
create_plan_ffth_(hipfftHandle * *plan, int *ISIGNp, int *Np, int *LOTp)
{
int ISIGN = *ISIGNp;
int N = *Np;
int LOT = *LOTp;

*plan = new hipfftHandle;
//hipfftHandle plan;

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Hip error: Failed to synchronize\n");
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
hipfftType hipfft_1 = HIPFFT_R2C;
hipfftType hipfft_2 = HIPFFT_C2R;
#else
hipfftType hipfft_1 = HIPFFT_D2Z;
hipfftType hipfft_2 = HIPFFT_Z2D;
#endif

embed[0] = 1;
stride   = LOT;
dist     = 1;

hipfftSafeCall(hipfftCreate(*plan));

//printf("CreatePlan hipfft\n","N=",N);
//printf("%s %d \n","plan=",plan);
//printf("%s %d \n","LOT=",LOT);
//printf("%s %d \n","ISIGN=",ISIGN);
//printf("%s %d \n","Np=",*Np);

if( ISIGN== -1 ){
  hipfftSafeCall(hipfftPlanMany(*plan, 1, &N,
                 embed, stride, dist, 
                 embed, stride, dist, 
                 hipfft_1, LOT));
  //hipfftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_D2Z, LOT));
}
else if( ISIGN== 1){
  hipfftSafeCall(hipfftPlanMany(*plan, 1, &N,
                 embed, stride, dist, 
                 embed, stride, dist, 
                 hipfft_2, LOT));
  //hipfftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_Z2D, LOT));
}
else {
  abort();
}

// // use our reusaable work area for the plan
// hipfftSetWorkArea(plan,planWorkspace); 

/*
if( ISIGN== -1 ){
  hipfftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_D2Z, LOT));
}
else if( ISIGN== 1){
  hipfftSafeCall(hipfftPlan1d(&plan, N, HIPFFT_Z2D, LOT));
}
else {
  abort();
}
*/

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Hip error: Failed to synchronize\n");
	return;	
}

//*PLANp=plan;
//fprintf(stderr, "create_plan_ffth_: plan-address = %p\n",*plan);

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

