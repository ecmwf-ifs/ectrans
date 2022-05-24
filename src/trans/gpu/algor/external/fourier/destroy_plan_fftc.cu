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

extern "C" void destroy_plan_fftc_(cufftHandle *PLANp) {
  cufftHandle plan = *PLANp;

  if (cudaDeviceSynchronize() != cudaSuccess) {
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }

  cufftSafeCall(cufftDestroy(plan));

  if (cudaDeviceSynchronize() != cudaSuccess) {
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }
}
