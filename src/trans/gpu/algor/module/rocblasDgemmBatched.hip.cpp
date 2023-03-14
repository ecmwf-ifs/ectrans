//
// Wrapper for hipblasDgemm function. 
//
// Alan Gray, NVIDIA
//
#include <iostream>
#include <stdio.h>
#include "hip/hip_runtime_api.h"
#include "rocblas.h"

using namespace std;

bool roc_alreadyAllocated_dgemm=false;
bool roc_alreadyAllocated_dgemm_handle=false;

double **d_Aarray_roc;
double **d_Barray_roc;
double **d_Carray_roc;

double **Aarray_roc;
double **Barray_roc;
double **Carray_roc;

//hipblasHandle_t handle_roc_dgemm;	
rocblas_handle handle_roc_dgemm;	

extern "C" void rocblasDgemmBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, int tda, const double *B, int ldb, int tdb, double beta, double *C, int ldc, int tdc, int batchCount)
{

  // printf("ROCBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
  //  hipblasStatus_t stat;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  //if (transa=='T' || transa=='t')	
  //  op_t1=HIPBLAS_OP_T;
  //if (transb=='T' || transb=='t')
  //  op_t2=HIPBLAS_OP_T;
  rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  if (transa=='T' || transa=='t') op_t1 = rocblas_operation_transpose;
  if (transb=='T' || transb=='t') op_t2 = rocblas_operation_transpose;

  //double **Aarray_roc = (double**) malloc(batchCount*sizeof(double*));
  //double **Barray_roc = (double**) malloc(batchCount*sizeof(double*));
  //double **Carray_roc = (double**) malloc(batchCount*sizeof(double*));

  if (!roc_alreadyAllocated_dgemm_handle){
    rocblas_status stat = rocblas_create_handle(&handle_roc_dgemm);
    if(stat == rocblas_status_success)
    {
      //cout << "status == rocblas_status_success" << endl;
    }
    else
    {
       cout << "rocblas failure: status = " << stat << endl;
    }
    //stat = hipblasCreate(&handle_roc_dgemm);
     //if (stat != HIPBLAS_STATUS_SUCCESS) {
     //   printf ("ROCBLAS initialization failed\n");
        //return EXIT_FAILURE;
    //}
  }
  roc_alreadyAllocated_dgemm_handle=true;

  if (!roc_alreadyAllocated_dgemm){
    hipError_t errcm1 = hipHostMalloc(&Aarray_roc,batchCount*sizeof(double*));
    hipError_t errcm2 = hipHostMalloc(&Barray_roc,batchCount*sizeof(double*));
    hipError_t errcm3 = hipHostMalloc(&Carray_roc,batchCount*sizeof(double*));
        
    hipError_t errcm4 = hipMalloc(&d_Aarray_roc,batchCount*sizeof(double*));
    hipError_t errcm5 = hipMalloc(&d_Barray_roc,batchCount*sizeof(double*));
    hipError_t errcm6 = hipMalloc(&d_Carray_roc,batchCount*sizeof(double*));
   }
  roc_alreadyAllocated_dgemm=true;

  int i;
  for(i=0;i<batchCount;i++){
    Aarray_roc[i]=(double*) &(A[i*lda*tda]);
    Barray_roc[i]=(double*) &(B[i*ldb*tdb]);
    Carray_roc[i]=(double*) &(C[i*ldc*tdc]);
  }

  hipError_t err1 = hipMemcpy(d_Aarray_roc,Aarray_roc,batchCount*sizeof(double*),hipMemcpyHostToDevice);
  hipError_t err2 = hipMemcpy(d_Barray_roc,Barray_roc,batchCount*sizeof(double*),hipMemcpyHostToDevice);
  hipError_t err3 = hipMemcpy(d_Carray_roc,Carray_roc,batchCount*sizeof(double*),hipMemcpyHostToDevice);

  //hipblasDgemmBatched(handle_roc_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray_roc,lda, (const double**) d_Barray_roc,ldb,&beta,(double**) d_Carray_roc,ldc,batchCount);
  rocblas_status status = rocblas_dgemm_batched(handle_roc_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray_roc,lda, (const double**) d_Barray_roc,ldb,&beta,(double**) d_Carray_roc,ldc,batchCount);
  if(status == rocblas_status_success)
    {
      //cout << "status == sgemm rocblas_status_success" << endl;
    }
    else
    {
       cout << "rocblas sgemm failure: status = " << status << endl;
    }
  // needed ?  
  //hipError_t err4 = hipMemcpy(Carray_roc,d_Carray_roc,batchCount*sizeof(double*),hipMemcpyDeviceToHost);
  //for(i=0;i<batchCount;i++){
  //  C[i]=*Carray_roc[i];
  //}
  hipDeviceSynchronize();
  
  //hipFree(Aarray_roc);
  //hipFree(Barray_roc);
  //hipFree(Carray_roc);
  
  //hipFree(d_Aarray_roc);
  //hipFree(d_Barray_roc);
  //hipFree(d_Carray_roc);
  //hipblasDestroy(handle_roc_dgemm);
  
}

extern "C" void rocblasDgemmStridedBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, long long tda, const double *B, int ldb, long long tdb, double beta, double *C, int ldc, long long tdc, int batchCount)
{


  // printf("ROCBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
  rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  if (transa=='T' || transa=='t') op_t1 = rocblas_operation_transpose;
  if (transb=='T' || transb=='t') op_t2 = rocblas_operation_transpose;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  //if (transa=='T' || transa=='t')       
  //  op_t1=HIPBLAS_OP_T;
  //if (transb=='T' || transb=='t')       
  //  op_t2=HIPBLAS_OP_T;

  if (!roc_alreadyAllocated_dgemm_handle){
    //hipblasCreate(&handle_roc_dgemm);
    rocblas_create_handle(&handle_roc_dgemm);
    roc_alreadyAllocated_dgemm_handle=true;
  }
  rocblas_dgemm_strided_batched(handle_roc_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,&beta,(double *) C,ldc,tdc,batchCount);
  //hipblasDgemmStridedBatched(handle_roc_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,&beta,(double *) C,ldc,tdc,batchCount);

}

extern "C" void rocblasDgemmBatched_finalize ()
{



  if (roc_alreadyAllocated_dgemm){
  
    hipFree(Aarray_roc);
    hipFree(Barray_roc);
    hipFree(Carray_roc);
    
    hipFree(d_Aarray_roc);
    hipFree(d_Barray_roc);
    hipFree(d_Carray_roc);
    if (roc_alreadyAllocated_dgemm_handle){
      //hipblasDestroy(handle_roc_dgemm);
      rocblas_destroy_handle(handle_roc_dgemm);
    }
    roc_alreadyAllocated_dgemm_handle=false;

  }
  roc_alreadyAllocated_dgemm=false;
  
}
