//
// Wrapper for hipblasSgemm function. 
//
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "hipblas.h" 


bool hip_alreadyAllocated_sgemm=false;
bool hip_alreadyAllocated_sgemm_handle=false;

float **d_Aarray_sgemm_hip;
float **d_Barray_sgemm_hip;
float **d_Carray_sgemm_hip;

float **Aarray_sgemm_hip;
float **Barray_sgemm_hip;
float **Carray_sgemm_hip;

hipblasHandle_t handle_hip_sgemm;	

extern "C" void hipblasSgemmBatched_wrapper (char transa, char transb, int m, int n,int k, float alpha, const float *A, int lda, int tda, const float *B, int ldb, int tdb, float beta, float *C, int ldc, int tdc, int batchCount)
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

  if (!hip_alreadyAllocated_sgemm_handle){
    hipblasCreate(&handle_hip_sgemm);
    hip_alreadyAllocated_sgemm_handle=true;
  }

  if (!hip_alreadyAllocated_sgemm){
    hipHostMalloc(&Aarray_sgemm_hip,batchCount*sizeof(float*),hipHostMallocNonCoherent);
    hipHostMalloc(&Barray_sgemm_hip,batchCount*sizeof(float*),hipHostMallocNonCoherent);
    hipHostMalloc(&Carray_sgemm_hip,batchCount*sizeof(float*),hipHostMallocNonCoherent);
 
    hipMalloc(&d_Aarray_sgemm_hip,batchCount*sizeof(float*));
    hipMalloc(&d_Barray_sgemm_hip,batchCount*sizeof(float*));
    hipMalloc(&d_Carray_sgemm_hip,batchCount*sizeof(float*));
    hip_alreadyAllocated_sgemm=true;
  }

  int i;
  for(i=0;i<batchCount;i++){
    Aarray_sgemm_hip[i]=(float*) &(A[i*lda*tda]);
    Barray_sgemm_hip[i]=(float*) &(B[i*ldb*tdb]);
    Carray_sgemm_hip[i]=(float*) &(C[i*ldc*tdc]);
  }
  hipMemcpy(d_Aarray_sgemm_hip,Aarray_sgemm_hip,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Barray_sgemm_hip,Barray_sgemm_hip,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Carray_sgemm_hip,Carray_sgemm_hip,batchCount*sizeof(float*),hipMemcpyHostToDevice);

  hipblasSgemmBatched(handle_hip_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float**) d_Aarray_sgemm_hip,lda, (const float**) d_Barray_sgemm_hip,ldb,&beta,(float**) d_Carray_sgemm_hip,ldc,batchCount);

  //printf("after sgemm\n");
  hipDeviceSynchronize();
  
  //hipFree(Aarray_sgemm_hip);
  //hipFree(Barray_sgemm_hip);
  //hipFree(Carray_sgemm_hip);
  
  //hipFree(d_Aarray_sgemm_hip);
  //hipFree(d_Barray_sgemm_hip);
  //hipFree(d_Carray_sgemm_hip);
  //hipblasDestroy(handle_hip_sgemm);
  
}

extern "C" void hipblasSgemmStridedBatched_wrapper (char transa, char transb, int m, int n,int k, float alpha, const float *A, int lda, long long tda, const float *B, int ldb, long long tdb, float beta, float *C, int ldc, long long tdc, int batchCount)
{


  printf("HIPBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
  exit;
 
  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')       
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')       
    op_t2=HIPBLAS_OP_T;

  if (!hip_alreadyAllocated_sgemm_handle){
    hipblasCreate(&handle_hip_sgemm);
    hip_alreadyAllocated_sgemm_handle=true;
  }
  hipblasSgemmStridedBatched(handle_hip_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float *) A,lda,tda, (const float *) B,ldb,tdb,&beta,(float*) C,ldc,tdc,batchCount);

}

extern "C" void hipblasSgemmBatched_finalize ()
{

  if (hip_alreadyAllocated_sgemm){
  
    hipFree(Aarray_sgemm_hip);
    hipFree(Barray_sgemm_hip);
    hipFree(Carray_sgemm_hip);
    
    hipFree(d_Aarray_sgemm_hip);
    hipFree(d_Barray_sgemm_hip);
    hipFree(d_Carray_sgemm_hip);

  }

  if (hip_alreadyAllocated_sgemm_handle){
    hipblasDestroy(handle_hip_sgemm);
  }
  
}
