//
// Wrapper for hipblasDgemm function. 
//
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "hipblas.h" 


bool alreadyAllocated_dgemm=false;
bool alreadyAllocated_dgemm_handle=false;

double **d_Aarray;
double **d_Barray;
double **d_Carray;

double **Aarray;
double **Barray;
double **Carray;

hipblasHandle_t handle_dgemm;	

extern "C" void cublasDgemmBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, int tda, const double *B, int ldb, int tdb, double beta, double *C, int ldc, int tdc, int batchCount)
{


  // printf("CUBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
    hipblasStatus_t stat;

 
  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')	
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=HIPBLAS_OP_T;


  //double **Aarray = (double**) malloc(batchCount*sizeof(double*));
  //double **Barray = (double**) malloc(batchCount*sizeof(double*));
  //double **Carray = (double**) malloc(batchCount*sizeof(double*));



  if (!alreadyAllocated_dgemm_handle){
     stat = hipblasCreate(&handle_dgemm);
     if (stat != HIPBLAS_STATUS_SUCCESS) {
        printf ("CUBLAS initialization failed\n");
        //return EXIT_FAILURE;
    }
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
  hipDeviceSynchronize();


  hipblasDgemmBatched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray,lda, (const double**) d_Barray,ldb,&beta,(double**) d_Carray,ldc,batchCount);

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

 
  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')       
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')       
    op_t2=HIPBLAS_OP_T;

  if (!alreadyAllocated_dgemm_handle){
    hipblasCreate(&handle_dgemm);
    alreadyAllocated_dgemm_handle=true;
  }
  hipblasDgemmStridedBatched(handle_dgemm,op_t1,op_t2,m,n,k,&alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,&beta,(double *) C,ldc,tdc,batchCount);

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
      hipblasDestroy(handle_dgemm);
    }
    alreadyAllocated_dgemm_handle=false;

  }
  alreadyAllocated_dgemm=false;
  
}
