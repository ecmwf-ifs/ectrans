#include "cufft.h"
#include "stdio.h"
extern "C" cufftDoubleComplex *create_storage_(int *Np) {
  int N = *Np;
  cufftDoubleComplex *data;
  /*cudaMalloc((void**)&data,sizeof(cufftDoubleComplex)*N);
  if (cudaGetLastError() != cudaSuccess){
       fprintf(stderr, "Cuda error: Failed to allocate\n");
       return 0;
  }
  return data;*/
  printf("%s %d \n", "sizeof(cufftDoubleComplex)=", sizeof(cufftDoubleComplex));
  printf("%s %d \n", "N=", N);
  if (cudaMalloc(&data, sizeof(cufftDoubleComplex) * N) == cudaSuccess) {
    printf("%s %X \n", "data ", data);
    return data;
  }
  fprintf(stderr, "Cuda error: Failed to allocate\n");
  return 0;
}

extern "C" void destroy_storage_(int *ptr) { cudaFree(ptr); }
