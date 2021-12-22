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

extern "C"
void
fft993_(cufftDoubleComplex *data, int *INCp, \
        int *JUMPp, int *Np, int *LOTp, int *ISIGNp)
{
cufftHandle plan;
int INC = *INCp;
int JUMP = *JUMPp;
int N = *Np;
int LOT = *LOTp;
int ISIGN = *ISIGNp;

/*
int RANK = 1;
int INEMBED[]={N};
int ISTRIDE = INC;
int IDIST = JUMP;
int ONEMBED[]={N/2+1};
int OSTRIDE = INC;
int ODIST = JUMP;

int NN[1] = {N};

printf("%s %d \n","sizeof(cufftDoubleComplex)=",sizeof(cufftDoubleComplex));
printf("%s %d \n","INC=",INC);
printf("%s %d \n","JUMP=",JUMP);
printf("%s %d \n","N=",N);
printf("%s %d \n","LOT=",LOT);
printf("%s %d \n","ISIGN=",ISIGN);
printf("%s %d \n","sizeof(cufftDoubleComplex)*(N/2+1)*LOT=",sizeof(cufftDoubleComplex)*(N/2+1)*LOT);
*/

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

if( ISIGN== -1 ){
  cufftSafeCall(cufftPlan1d(&plan, N, CUFFT_D2Z, LOT));
  /*
  cufftSafeCall(cufftPlanMany(&plan, RANK, NN, INEMBED, ISTRIDE, IDIST, ONEMBED, OSTRIDE, ODIST, CUFFT_D2Z, LOT ));
  */
  /* Use the CUFFT plan to transform the signal in place. */
  cufftSafeCall(cufftExecD2Z(plan, (cufftDoubleReal*)data, data));
}
else if( ISIGN== 1){
  cufftSafeCall(cufftPlan1d(&plan, N, CUFFT_Z2D, LOT));
  /*
  cufftSafeCall(cufftPlanMany(&plan, RANK, NN, ONEMBED, OSTRIDE, ODIST, INEMBED, ISTRIDE, IDIST, CUFFT_Z2D, LOT ));
  */
  /* Use the CUFFT plan to transform the signal in place. */
  cufftSafeCall(cufftExecZ2D(plan, data, (cufftDoubleReal*)data));
}
else {
  abort();
}


if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

cufftDestroy(plan);

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

}

