#include "hipfft.h"
#include "stdio.h"
extern "C"
hipfftDoubleComplex *create_storage_(int *Np)
{
   int N = *Np;
   hipfftDoubleComplex *data;
   /*hipMalloc((void**)&data,sizeof(hipfftDoubleComplex)*N);
   if (hipGetLastError() != hipSuccess){
        fprintf(stderr, "Cuda error: Failed to allocate\n");
        return 0;
   } 
   return data;*/
   printf("%s %d \n","sizeof(hipfftDoubleComplex)=",sizeof(hipfftDoubleComplex));
   printf("%s %d \n","N=",N);
   if (hipMalloc(&data, sizeof(hipfftDoubleComplex)*N) == hipSuccess){ 
        printf("%s %X \n","data ",data);
        return data;
   }
   fprintf(stderr, "Cuda error: Failed to allocate\n");
   return 0;
}

extern "C"
void destroy_storage_(int *ptr)
{
   hipFree(ptr);
}
