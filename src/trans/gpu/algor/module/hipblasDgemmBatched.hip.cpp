//
// Wrapper for hipblasDgemm function. 
//
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "hipblas.h" 


bool hip_alreadyAllocated_dgemm=false;
bool hip_alreadyAllocated_dgemm_handle=false;

double **d_Aarray_hip;
double **d_Barray_hip;
double **d_Carray_hip;

double **Aarray_hip;
double **Barray_hip;
double **Carray_hip;

hipblasHandle_t handle_hip_dgemm;	

extern "C" void hipblasDgemmBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, int tda, const double *B, int ldb, int tdb, double beta, double *C, int ldc, int tdc, int batchCount)
{


  // printf("HIPBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
    hipblasStatus_t stat;

 
  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')	
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=HIPBLAS_OP_T;


  //double **Aarray_hip = (double**) malloc(batchCount*sizeof(double*));
  //double **Barray_hip = (double**) malloc(batchCount*sizeof(double*));
  //double **Carray_hip = (double**) malloc(batchCount*sizeof(double*));



  if (!hip_alreadyAllocated_dgemm_handle){
     stat = hipblasCreate(&handle_hip_dgemm);
     if (stat != HIPBLAS_STATUS_SUCCESS) {
        printf ("HIPBLAS initialization failed\n");
        //return EXIT_FAILURE;
    }
  }
  hip_alreadyAllocated_dgemm_handle=true;

  if (!hip_alreadyAllocated_dgemm){
    hipError_t errcm1 = hipHostMalloc(&Aarray_hip,batchCount*sizeof(double*),hipHostMallocNonCoherent);
    hipError_t errcm2 = hipHostMalloc(&Barray_hip,batchCount*sizeof(double*),hipHostMallocNonCoherent);
    hipError_t errcm3 = hipHostMalloc(&Carray_hip,batchCount*sizeof(double*),hipHostMallocNonCoherent);
        
    hipError_t errcm4 = hipMalloc(&d_Aarray_hip,batchCount*sizeof(double*));
    hipError_t errcm5 = hipMalloc(&d_Barray_hip,batchCount*sizeof(double*));
    hipError_t errcm6 = hipMalloc(&d_Carray_hip,batchCount*sizeof(double*));
   }
  hip_alreadyAllocated_dgemm=true;

  int i;
  for(i=0;i<batchCount;i++){
    Aarray_hip[i]=(double*) &(A[i*lda*tda]);
    Barray_hip[i]=(double*) &(B[i*ldb*tdb]);
    Carray_hip[i]=(double*) &(C[i*ldc*tdc]);
  }

  hipError_t err1 = hipMemcpy(d_Aarray_hip,Aarray_hip,batchCount*sizeof(double*),hipMemcpyHostToDevice);
  hipError_t err2 = hipMemcpy(d_Barray_hip,Barray_hip,batchCount*sizeof(double*),hipMemcpyHostToDevice);
  hipError_t err3 = hipMemcpy(d_Carray_hip,Carray_hip,batchCount*sizeof(double*),hipMemcpyHostToDevice);
  hipDeviceSynchronize();


  hipblasDgemmBatched(handle_hip_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray_hip,lda, (const double**) d_Barray_hip,ldb,&beta,(double**) d_Carray_hip,ldc,batchCount);

  hipDeviceSynchronize();
  
  //hipFree(Aarray_hip);
  //hipFree(Barray_hip);
  //hipFree(Carray_hip);
  
  //hipFree(d_Aarray_hip);
  //hipFree(d_Barray_hip);
  //hipFree(d_Carray_hip);
  //hipblasDestroy(handle_hip_dgemm);
  
  
}

extern "C" void hipblasDgemmStridedBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, long long tda, const double *B, int ldb, long long tdb, double beta, double *C, int ldc, long long tdc, int batchCount)
{


  // printf("HIPBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);

 
  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')       
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')       
    op_t2=HIPBLAS_OP_T;

  if (!hip_alreadyAllocated_dgemm_handle){
    hipblasCreate(&handle_hip_dgemm);
    hip_alreadyAllocated_dgemm_handle=true;
  }
  hipblasDgemmStridedBatched(handle_hip_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,&beta,(double *) C,ldc,tdc,batchCount);

}

extern "C" void hipblasDgemmBatched_finalize ()
{



  if (hip_alreadyAllocated_dgemm){
  
    hipFree(Aarray_hip);
    hipFree(Barray_hip);
    hipFree(Carray_hip);
    
    hipFree(d_Aarray_hip);
    hipFree(d_Barray_hip);
    hipFree(d_Carray_hip);
    if (hip_alreadyAllocated_dgemm_handle){
      hipblasDestroy(handle_hip_dgemm);
    }
    hip_alreadyAllocated_dgemm_handle=false;

  }
  hip_alreadyAllocated_dgemm=false;
  
}
