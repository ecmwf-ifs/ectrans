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

bool alreadyAllocated_dgemm=false;
bool alreadyAllocated_dgemm_handle=false;

double **d_Aarray;
double **d_Barray;
double **d_Carray;

double **Aarray;
double **Barray;
double **Carray;

//hipblasHandle_t handle_dgemm;	
rocblas_handle handle_dgemm;	

extern "C" void cublasDgemmBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, int tda, const double *B, int ldb, int tdb, double beta, double *C, int ldc, int tdc, int batchCount)
{

  // printf("CUBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
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

  //double **Aarray = (double**) malloc(batchCount*sizeof(double*));
  //double **Barray = (double**) malloc(batchCount*sizeof(double*));
  //double **Carray = (double**) malloc(batchCount*sizeof(double*));

  if (!alreadyAllocated_dgemm_handle){
    rocblas_status stat = rocblas_create_handle(&handle_dgemm);
    if(stat == rocblas_status_success)
    {
      //cout << "status == rocblas_status_success" << endl;
    }
    else
    {
       cout << "rocblas failure: status = " << stat << endl;
    }
    //stat = hipblasCreate(&handle_dgemm);
     //if (stat != HIPBLAS_STATUS_SUCCESS) {
     //   printf ("CUBLAS initialization failed\n");
        //return EXIT_FAILURE;
    //}
  }
  alreadyAllocated_dgemm_handle=true;

  if (!alreadyAllocated_dgemm){
    hipError_t errcm1 = hipHostMalloc(&Aarray,batchCount*sizeof(double*));
    hipError_t errcm2 = hipHostMalloc(&Barray,batchCount*sizeof(double*));
    hipError_t errcm3 = hipHostMalloc(&Carray,batchCount*sizeof(double*));
        
    hipError_t errcm4 = hipMalloc(&d_Aarray,batchCount*sizeof(double*));
    hipError_t errcm5 = hipMalloc(&d_Barray,batchCount*sizeof(double*));
    hipError_t errcm6 = hipMalloc(&d_Carray,batchCount*sizeof(double*));
   }
  alreadyAllocated_dgemm=true;

  int i;
  for(i=0;i<batchCount;i++){
    Aarray[i]=(double*) &(A[i*lda*tda]);
    Barray[i]=(double*) &(B[i*ldb*tdb]);
    Carray[i]=(double*) &(C[i*ldc*tdc]);
  }

  hipError_t err1 = hipMemcpy(d_Aarray,Aarray,batchCount*sizeof(double*),hipMemcpyHostToDevice);
  hipError_t err2 = hipMemcpy(d_Barray,Barray,batchCount*sizeof(double*),hipMemcpyHostToDevice);
  hipError_t err3 = hipMemcpy(d_Carray,Carray,batchCount*sizeof(double*),hipMemcpyHostToDevice);

  //hipblasDgemmBatched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray,lda, (const double**) d_Barray,ldb,&beta,(double**) d_Carray,ldc,batchCount);
  rocblas_status status = rocblas_dgemm_batched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray,lda, (const double**) d_Barray,ldb,&beta,(double**) d_Carray,ldc,batchCount);
  if(status == rocblas_status_success)
    {
      //cout << "status == sgemm rocblas_status_success" << endl;
    }
    else
    {
       cout << "rocblas sgemm failure: status = " << status << endl;
    }
  // needed ?  
  //hipError_t err4 = hipMemcpy(Carray,d_Carray,batchCount*sizeof(double*),hipMemcpyDeviceToHost);
  //for(i=0;i<batchCount;i++){
  //  C[i]=*Carray[i];
  //}
  hipDeviceSynchronize();
  
  //hipFree(Aarray);
  //hipFree(Barray);
  //hipFree(Carray);
  
  //hipFree(d_Aarray);
  //hipFree(d_Barray);
  //hipFree(d_Carray);
  //hipblasDestroy(handle_dgemm);
  
}

extern "C" void cublasDgemmStridedBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, long long tda, const double *B, int ldb, long long tdb, double beta, double *C, int ldc, long long tdc, int batchCount)
{


  // printf("CUBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
  rocblas_operation op_t1 = rocblas_operation_none, op_t2 = rocblas_operation_none;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  if (transa=='T' || transa=='t') op_t1 = rocblas_operation_transpose;
  if (transb=='T' || transb=='t') op_t2 = rocblas_operation_transpose;
  //hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;
  //if (transa=='T' || transa=='t')       
  //  op_t1=HIPBLAS_OP_T;
  //if (transb=='T' || transb=='t')       
  //  op_t2=HIPBLAS_OP_T;

  if (!alreadyAllocated_dgemm_handle){
    //hipblasCreate(&handle_dgemm);
    rocblas_create_handle(&handle_dgemm);
    alreadyAllocated_dgemm_handle=true;
  }
  rocblas_dgemm_strided_batched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,&beta,(double *) C,ldc,tdc,batchCount);
  //hipblasDgemmStridedBatched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,&beta,(double *) C,ldc,tdc,batchCount);

}

extern "C" void cublasDgemmBatched_finalize ()
{



  if (alreadyAllocated_dgemm){
  
    hipFree(Aarray);
    hipFree(Barray);
    hipFree(Carray);
    
    hipFree(d_Aarray);
    hipFree(d_Barray);
    hipFree(d_Carray);
    if (alreadyAllocated_dgemm_handle){
      //hipblasDestroy(handle_dgemm);
      rocblas_destroy_handle(handle_dgemm);
    }
    alreadyAllocated_dgemm_handle=false;

  }
  alreadyAllocated_dgemm=false;
  
}
