#define cufftSafeCall(err) __cufftSafeCall(err, __FILE__, __LINE__)
#include "cufft.h"
#include "stdio.h"
static const char *_cudaGetErrorEnum(cufftResult error) {
  switch (error) {
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

inline void __cufftSafeCall(cufftResult err, const char *file, const int line) {
  if (CUFFT_SUCCESS != err) {
    fprintf(stderr, "CUFFT error at 1\n");
    fprintf(stderr, "CUFFT error in file '%s'\n", __FILE__);
    fprintf(stderr, "CUFFT error at 2\n");
    /*fprintf(stderr, "CUFFT error line '%s'\n",__LINE__);*/
    fprintf(stderr, "CUFFT error at 3\n");
    /*fprintf(stderr, "CUFFT error in file '%s', line %d\n %s\nerror %d:
    %s\nterminating!\n",__FILE__, __LINE__,err, \
    _cudaGetErrorEnum(err)); \*/
    fprintf(stderr, "CUFFT error %d: %s\nterminating!\n", err,
            _cudaGetErrorEnum(err));
    cudaDeviceReset();
    return;
  }
}

void *planWorkspace;
static int currentWorkspaceSize = 0;

extern "C" void create_plan_fftc_(cufftHandle *PLANp, int *ISIGNp, int *Np,
                                  int *LOTp, int *stridep, int *plan_size) {
  int ISIGN = *ISIGNp;
  int N = *Np;
  int LOT = *LOTp;
  int stride = *stridep;

  cufftHandle plan;

  if (cudaDeviceSynchronize() != cudaSuccess) {
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }

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
  dist = 1;

  cufftSafeCall(cufftCreate(&plan));

  // Disable auto allocation
  cufftSetAutoAllocation(plan, false);

  // printf("CreatePlan cuFFT\n","N=",N);
  // printf("%s %d \n","plan=",plan);
  // printf("%s %d \n","LOT=",LOT);
  // printf("%s %d \n","ISIGN=",ISIGN);
  // printf("%s %d \n","Np=",*Np);

  if (ISIGN == -1) {
    cufftSafeCall(cufftPlanMany(&plan, 1, &N, embed, stride, dist, embed,
                                stride, dist, cufft_1, LOT));
  } else if (ISIGN == 1) {
    cufftSafeCall(cufftPlanMany(&plan, 1, &N, embed, stride, dist, embed,
                                stride, dist, cufft_2, LOT));
  } else {
    abort();
  }

  // get size used by this plan
  size_t thisWorkplanSize;
  cufftGetSize(plan, &thisWorkplanSize);

  // check if this the work space is sufficiently large
  if (thisWorkplanSize > currentWorkspaceSize) {
    cudaDeviceSynchronize();
    cudaFree(planWorkspace);
    cudaMalloc(&planWorkspace, thisWorkplanSize);
    currentWorkspaceSize = thisWorkplanSize;
  }

  if (cudaDeviceSynchronize() != cudaSuccess) {
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }

  *PLANp = plan;
  *plan_size = thisWorkplanSize;

  return;
}
