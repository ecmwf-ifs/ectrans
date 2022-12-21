//
// Wrapper for hipblasSgemm function. 
//
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "hipblas.h" 


bool alreadyAllocated_sgemm=false;
bool alreadyAllocated_sgemm_handle=false;

float **d_Aarray_sgemm;
float **d_Barray_sgemm;
float **d_Carray_sgemm;

float **Aarray_sgemm;
float **Barray_sgemm;
float **Carray_sgemm;

hipblasHandle_t handle_sgemm;	

extern "C" void cublasSgemmBatched_wrapper (char transa, char transb, int m, int n,int k, float alpha, const float *A, int lda, int tda, const float *B, int ldb, int tdb, float beta, float *C, int ldc, int tdc, int batchCount)
{

   printf("HIPBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
   //exit;
  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')		
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')		
    op_t2=HIPBLAS_OP_T;

  //float **Aarray_sgemm = (float**) malloc(batchCount*sizeof(float*));
  //float **Barray_sgemm = (float**) malloc(batchCount*sizeof(float*));
  //float **Carray_sgemm = (float**) malloc(batchCount*sizeof(float*));

  if (!alreadyAllocated_sgemm_handle){
    hipblasCreate(&handle_sgemm);
    alreadyAllocated_sgemm_handle=true;
  }

  if (!alreadyAllocated_sgemm){
    hipHostMalloc(&Aarray_sgemm,batchCount*sizeof(float*));
    hipHostMalloc(&Barray_sgemm,batchCount*sizeof(float*));
    hipHostMalloc(&Carray_sgemm,batchCount*sizeof(float*));
    alreadyAllocated_sgemm=true;
  }

  hipMalloc(&d_Aarray_sgemm,batchCount*sizeof(float*));
  hipMalloc(&d_Barray_sgemm,batchCount*sizeof(float*));
  hipMalloc(&d_Carray_sgemm,batchCount*sizeof(float*));

  int i;
  for(i=0;i<batchCount;i++){
    Aarray_sgemm[i]=(float*) &(A[i*lda*tda]);
    Barray_sgemm[i]=(float*) &(B[i*ldb*tdb]);
    Carray_sgemm[i]=(float*) &(C[i*ldc*tdc]);
  }
  hipMemcpy(d_Aarray_sgemm,Aarray_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Barray_sgemm,Barray_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Carray_sgemm,Carray_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);

  hipblasSgemmBatched(handle_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float**) d_Aarray_sgemm,lda, (const float**) d_Barray_sgemm,ldb,&beta,(float**) d_Carray_sgemm,ldc,batchCount);

  //printf("after sgemm\n");
  hipDeviceSynchronize();
  
  //hipFree(Aarray_sgemm);
  //hipFree(Barray_sgemm);
  //hipFree(Carray_sgemm);
  
  hipFree(d_Aarray_sgemm);
  hipFree(d_Barray_sgemm);
  hipFree(d_Carray_sgemm);
  //hipblasDestroy(handle_sgemm);
  
}

extern "C" void cublasSgemmStridedBatched_wrapper (char transa, char transb, int m, int n,int k, float alpha, const float *A, int lda, long long tda, const float *B, int ldb, long long tdb, float beta, float *C, int ldc, long long tdc, int batchCount)
{


  printf("HIPBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
  exit;
 
  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')       
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')       
    op_t2=HIPBLAS_OP_T;

  if (!alreadyAllocated_sgemm_handle){
    hipblasCreate(&handle_sgemm);
    alreadyAllocated_sgemm_handle=true;
  }
  hipblasSgemmStridedBatched(handle_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float *) A,lda,tda, (const float *) B,ldb,tdb,&beta,(float*) C,ldc,tdc,batchCount);

}

extern "C" void cublasSgemmBatched_finalize ()
{

  if (alreadyAllocated_sgemm){
  
    hipFree(Aarray_sgemm);
    hipFree(Barray_sgemm);
    hipFree(Carray_sgemm);
    
    hipFree(d_Aarray_sgemm);
    hipFree(d_Barray_sgemm);
    hipFree(d_Carray_sgemm);

  }

  if (alreadyAllocated_sgemm_handle){
    hipblasDestroy(handle_sgemm);
  }
  
}
