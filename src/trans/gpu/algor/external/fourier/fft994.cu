#include "cufft.h"
#include "stdio.h"
extern "C"
void
fft994_(cufftDoubleComplex *data_h, int *INCp, \
        int *JUMPp, int *Np, int *LOTp, int *ISIGNp)
{
cufftHandle plan;
cufftDoubleComplex *data;
int INC = *INCp;
int JUMP = *JUMPp;
int N = *Np;
int LOT = *LOTp;
int ISIGN = *ISIGNp;
int RANK = 1;

int INEMBED[]={N};
int ISTRIDE = INC;
int IDIST = JUMP;
int ONEMBED[]={N/2+1};
int OSTRIDE = INC;
int ODIST = JUMP;

int NN[1] = {N};
/*
printf("%s %d \n","sizeof(cufftDoubleComplex)=",sizeof(cufftDoubleComplex));
printf("%s %d \n","INC=",INC);
printf("%s %d \n","JUMP=",JUMP);
printf("%s %d \n","N=",N);
printf("%s %d \n","LOT=",LOT);
printf("%s %d \n","ISIGN=",ISIGN);
printf("%s %d \n","sizeof(cufftDoubleComplex)*(N/2+1)*LOT=",sizeof(cufftDoubleComplex)*(N/2+1)*LOT);
*/

cudaMalloc((void**)&data, sizeof(cufftDoubleComplex)*(N/2+1)*LOT);
if (cudaGetLastError() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to allocate\n");
	return;	
}

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

cudaMemcpy( data, data_h, sizeof(cufftDoubleComplex)*(N/2+1)*LOT, cudaMemcpyHostToDevice );

if( ISIGN== -1 ){
  /*
  if (cufftPlan1d(&plan, N, CUFFT_D2Z, LOT) != CUFFT_SUCCESS){
	  fprintf(stderr, "CUFFT error(DIR): Plan creation failed");
	  return;	
  }	
  */
  if(cufftPlanMany(&plan, RANK, NN, INEMBED, \
        ISTRIDE, IDIST, ONEMBED, OSTRIDE, \
        ODIST, CUFFT_D2Z, LOT ) != CUFFT_SUCCESS){
        fprintf(stderr, "CUFFT error(DIR): Plan creation failed");
        return;
  }
  /* Use the CUFFT plan to transform the signal in place. */
  if (cufftExecD2Z(plan, (cufftDoubleReal*)data, data) != CUFFT_SUCCESS){
	fprintf(stderr, "CUFFT error(DIR): ExecD2Z failed");
	return;	
  }
}
else if( ISIGN== 1){
  /*
  if (cufftPlan1d(&plan, N, CUFFT_Z2D, LOT) != CUFFT_SUCCESS){
	  fprintf(stderr, "CUFFT error(INV): Plan creation failed");
	  return;	
  }	
  */
  if(cufftPlanMany(&plan, RANK, NN, ONEMBED, \
        OSTRIDE, ODIST, INEMBED, ISTRIDE, \
        IDIST, CUFFT_Z2D, LOT ) != CUFFT_SUCCESS){
        fprintf(stderr, "CUFFT error(INV): Plan creation failed");
        return;
  }
  /* Use the CUFFT plan to transform the signal in place. */
  if (cufftExecZ2D(plan, data, (cufftDoubleReal*)data) != CUFFT_SUCCESS){
	fprintf(stderr, "CUFFT error(INV): ExecZ2D failed");
	return;	
  }
}
else {
  abort();
}


if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

cudaMemcpy( data_h, data, sizeof(cufftDoubleComplex)*(N/2+1)*LOT, cudaMemcpyDeviceToHost );
cufftDestroy(plan);
cudaFree(data);

if (cudaDeviceSynchronize() != cudaSuccess){
	fprintf(stderr, "Cuda error: Failed to synchronize\n");
	return;	
}

}
