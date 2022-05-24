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

extern void *planWorkspace;

extern "C" void
#ifdef TRANS_SINGLE
execute_plan_fftc_(cufftHandle *PLANp, int *ISIGNp, cufftComplex *data_in,
                   cufftComplex *data_out)
#else
execute_plan_fftc_(cufftHandle *PLANp, int *ISIGNp, cufftDoubleComplex *data_in,
                   cufftDoubleComplex *data_out)
#endif
{
  cufftHandle plan = *PLANp;
  int ISIGN = *ISIGNp;

  cufftSafeCall(cufftSetWorkArea(plan, planWorkspace));

  if (ISIGN == -1) {
#ifdef TRANS_SINGLE
    cufftSafeCall(cufftExecR2C(plan, (cufftReal *)data_in, data_out));
#else
    cufftSafeCall(cufftExecD2Z(plan, (cufftDoubleReal *)data_in, data_out));
#endif
  } else if (ISIGN == 1) {
#ifdef TRANS_SINGLE
    cufftSafeCall(cufftExecC2R(plan, data_in, (cufftReal *)data_out));
#else
    cufftSafeCall(cufftExecZ2D(plan, data_in, (cufftDoubleReal *)data_out));
#endif
  } else {
    abort();
  }
}
