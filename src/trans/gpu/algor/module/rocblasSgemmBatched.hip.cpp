//
// Wrapper for hipblasSgemm function. 
//
#include <iostream>
#include <stdio.h>
#include "hip/hip_runtime_api.h"
#include "rocblas.h"
//#include "hipblas.h" 

using namespace std;

bool roc_alreadyAllocated_sgemm=false;
bool roc_alreadyAllocated_sgemm_handle=false;

float **d_Aarray_roc_sgemm;
float **d_Barray_roc_sgemm;
float **d_Carray_roc_sgemm;

float **Aarray_roc_sgemm;
float **Barray_roc_sgemm;
float **Carray_roc_sgemm;

//hipblasHandle_t handle_roc_sgemm;	
rocblas_handle handle_roc_sgemm;

extern "C" void rocblasSgemmBatched_wrapper (char transa, char transb, int m, int n,int k, float alpha, const float *A, int lda, int tda, const float *B, int ldb, int tdb, float beta, float *C, int ldc, int tdc, int batchCount)
{

   //TCo79 default	ROCBLAS m=4,n=80,k=40,batchcount=80
   //printf("ROCBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
   //exit;
  rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  if (transa=='T' || transa=='t') 
    op_t1 = rocblas_operation_transpose;
    //op_t1=HIPBLAS_OP_T;
  if (transb=='T' || transb=='t') 
    op_t2 = rocblas_operation_transpose;
    //op_t2=HIPBLAS_OP_T;

  //float **Aarray_roc_sgemm = (float**) malloc(batchCount*sizeof(float*));
  //float **Barray_roc_sgemm = (float**) malloc(batchCount*sizeof(float*));
  //float **Carray_roc_sgemm = (float**) malloc(batchCount*sizeof(float*));

  if (!roc_alreadyAllocated_sgemm_handle){
    //hipblasCreate(&handle_roc_sgemm);
    rocblas_status stat = rocblas_create_handle(&handle_roc_sgemm);
    if(stat == rocblas_status_success)
    { 
       //cout << "status create == rocblas_status_success" << endl;
    }
    else
    {
       cout << "rocblas create failure: status = " << stat << endl;
    }
  }
  roc_alreadyAllocated_sgemm_handle=true;

  if (!roc_alreadyAllocated_sgemm){
    hipHostMalloc(&Aarray_roc_sgemm,batchCount*sizeof(float*));
    hipHostMalloc(&Barray_roc_sgemm,batchCount*sizeof(float*));
    hipHostMalloc(&Carray_roc_sgemm,batchCount*sizeof(float*));
  }
  roc_alreadyAllocated_sgemm=true;

  hipMalloc(&d_Aarray_roc_sgemm,batchCount*sizeof(float*));
  hipMalloc(&d_Barray_roc_sgemm,batchCount*sizeof(float*));
  hipMalloc(&d_Carray_roc_sgemm,batchCount*sizeof(float*));

  int i;
  for(i=0;i<batchCount;i++){
    Aarray_roc_sgemm[i]=(float*) &(A[i*lda*tda]);
    Barray_roc_sgemm[i]=(float*) &(B[i*ldb*tdb]);
    Carray_roc_sgemm[i]=(float*) &(C[i*ldc*tdc]);
  }
  hipMemcpy(d_Aarray_roc_sgemm,Aarray_roc_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Barray_roc_sgemm,Barray_roc_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);
  hipMemcpy(d_Carray_roc_sgemm,Carray_roc_sgemm,batchCount*sizeof(float*),hipMemcpyHostToDevice);

  //for(i=0;i<batchCount;i++){
  //   cout << "before first element matrix for batch " << i << "val=" << *Carray_roc_sgemm[i] << endl;
  //}

  //hipblasSgemmBatched(handle_roc_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float**) d_Aarray_roc_sgemm,lda, (const float**) d_Barray_roc_sgemm,ldb,&beta,(float**) d_Carray_roc_sgemm,ldc,batchCount);
  rocblas_status status = rocblas_sgemm_batched(handle_roc_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float**) d_Aarray_roc_sgemm,lda, (const float**) d_Barray_roc_sgemm,ldb,&beta,(float**) d_Carray_roc_sgemm,ldc,batchCount);
  if(status == rocblas_status_success)
   {
      //cout << "status sgemm == rocblas_status_success" << endl;
   }
   else
   {
      cout << "rocblas sgemm failure: status = " << status << endl;
   }
  // need this for debugging
  // d_Carray_roc_sgemm maps to the memory space on the device occupied by the corresponding Fortran device array_roc
  //hipMemcpy(Carray_roc_sgemm,d_Carray_roc_sgemm,batchCount*sizeof(float*),hipMemcpyDeviceToHost);
  //for(i=0;i<batchCount;i++){
  //   cout << "after first element matrix for batch " << i << "val=" << *Carray_roc_sgemm[i] << endl;
  //}
  hipDeviceSynchronize();
  //printf("after sgemm\n");
  
  //hipFree(Aarray_roc_sgemm);
  //hipFree(Barray_roc_sgemm);
  //hipFree(Carray_roc_sgemm);
  
  hipFree(d_Aarray_roc_sgemm);
  hipFree(d_Barray_roc_sgemm);
  hipFree(d_Carray_roc_sgemm);
  //hipblasDestroy(handle_roc_sgemm);
  
}

extern "C" void rocblasSgemmStridedBatched_wrapper (char transa, char transb, int m, int n,int k, float alpha, const float *A, int lda, long long tda, const float *B, int ldb, long long tdb, float beta, float *C, int ldc, long long tdc, int batchCount)
{

  //printf("ROCBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
  rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  if (transa=='T' || transa=='t') op_t1 = rocblas_operation_transpose;
  //  op_t1=HIPBLAS_OP_T;
  if (transb=='T' || transb=='t') op_t2 = rocblas_operation_transpose;
  //  op_t2=HIPBLAS_OP_T;

  if (!roc_alreadyAllocated_sgemm_handle){
    rocblas_create_handle(&handle_roc_sgemm);
    //hipblasCreate(&handle_roc_sgemm);
    roc_alreadyAllocated_sgemm_handle=true;
  }
  //hipblasSgemmStridedBatched(handle_roc_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float *) A,lda,tda, (const float *) B,ldb,tdb,&beta,(float*) C,ldc,tdc,batchCount);
  rocblas_sgemm_strided_batched(handle_roc_sgemm,op_t1,op_t2,m,n,k,&alpha,(const float *) A,lda,tda, (const float *) B,ldb,tdb,&beta,(float*) C,ldc,tdc,batchCount);

}

extern "C" void rocblasSgemmBatched_finalize ()
{

  if (roc_alreadyAllocated_sgemm){
  
    hipFree(Aarray_roc_sgemm);
    hipFree(Barray_roc_sgemm);
    hipFree(Carray_roc_sgemm);
    
    hipFree(d_Aarray_roc_sgemm);
    hipFree(d_Barray_roc_sgemm);
    hipFree(d_Carray_roc_sgemm);

  }

  if (roc_alreadyAllocated_sgemm_handle){
    rocblas_destroy_handle(handle_roc_sgemm);
    //hipblasDestroy(handle_roc_sgemm);
  }
  
}
