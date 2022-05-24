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

#define CUFFT_CHECK(e) { \
	cufftResult_t err = (e); \
	if (err != CUFFT_SUCCESS) \
	{ \
		fprintf(stderr, "CUFFT error: %s, line %d, %s: %s\n", \
			__FILE__, __LINE__, #e, _cudaGetErrorEnum(err)); \
		exit(EXIT_FAILURE); \
	} \
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

  CUFFT_CHECK(cufftSetWorkArea(plan, planWorkspace));

  if (ISIGN == -1) {
#ifdef TRANS_SINGLE
    CUFFT_CHECK(cufftExecR2C(plan, (cufftReal *)data_in, data_out));
#else
    CUFFT_CHECK(cufftExecD2Z(plan, (cufftDoubleReal *)data_in, data_out));
#endif
  } else if (ISIGN == 1) {
#ifdef TRANS_SINGLE
    CUFFT_CHECK(cufftExecC2R(plan, data_in, (cufftReal *)data_out));
#else
    CUFFT_CHECK(cufftExecZ2D(plan, data_in, (cufftDoubleReal *)data_out));
#endif
  } else {
    abort();
  }
}
