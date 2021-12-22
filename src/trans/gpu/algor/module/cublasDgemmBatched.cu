//
// Wrapper for cublasDgemm function. 
//
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "cublas_v2.h" 


bool alreadyAllocated=false;

double **d_Aarray;
double **d_Barray;
double **d_Carray;

double **Aarray;
double **Barray;
double **Carray;

cublasHandle_t handle;	

extern "C" void cublasDgemmBatched_wrapper (char transa, char transb, int m, int n,int k, double alpha, const double *A, int lda, int tda, const double *B, int ldb, int tdb, double beta, double *C, int ldc, int tdc, int batchCount)
{


  // printf("CUBLAS m=%d,n=%d,k=%d,batchcount=%d\n",m,n,k,batchCount);
    cublasStatus_t stat;

 
  cublasOperation_t op_t1=CUBLAS_OP_N, op_t2=CUBLAS_OP_N;

  if (transa=='T' || transa=='t')	
    op_t1=CUBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=CUBLAS_OP_T;


  //double **Aarray = (double**) malloc(batchCount*sizeof(double*));
  //double **Barray = (double**) malloc(batchCount*sizeof(double*));
  //double **Carray = (double**) malloc(batchCount*sizeof(double*));


  if (!alreadyAllocated){

     stat = cublasCreate(&handle);
     if (stat != CUBLAS_STATUS_SUCCESS) {
        printf ("CUBLAS initialization failed\n");
        //return EXIT_FAILURE;
    }
    printf("cublascreate return code : %d\n",stat);

    cudaError_t errcm1 = cudaMallocHost(&Aarray,batchCount*sizeof(double*));
    cudaError_t errcm2 = cudaMallocHost(&Barray,batchCount*sizeof(double*));
    cudaError_t errcm3 = cudaMallocHost(&Carray,batchCount*sizeof(double*));
        
    cudaError_t errcm4 = cudaMalloc(&d_Aarray,batchCount*sizeof(double*));
    cudaError_t errcm5 = cudaMalloc(&d_Barray,batchCount*sizeof(double*));
    cudaError_t errcm6 = cudaMalloc(&d_Carray,batchCount*sizeof(double*));
 
    printf("switched alreadyAllocated to true\n");
    printf("Allocation statuses : %d %d %d %d %d %d\n", errcm1, errcm2, errcm3, errcm4, errcm5, errcm6 );
    alreadyAllocated=true;
  }

  int i;
  for(i=0;i<batchCount;i++){
    Aarray[i]=(double*) &(A[i*lda*tda]);
    Barray[i]=(double*) &(B[i*ldb*tdb]);
    Carray[i]=(double*) &(C[i*ldc*tdc]);
  }

  cudaError_t err1 = cudaMemcpy(d_Aarray,Aarray,batchCount*sizeof(double*),cudaMemcpyHostToDevice);
  cudaError_t err2 = cudaMemcpy(d_Barray,Barray,batchCount*sizeof(double*),cudaMemcpyHostToDevice);
  cudaError_t err3 = cudaMemcpy(d_Carray,Carray,batchCount*sizeof(double*),cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();

  printf("made it to the call to DgemmBatched ... are we already allocated? %d and err codes : %d %d %d \n",alreadyAllocated, err1, err2, err3);
  printf("batchCount etc : %d \n%d %d %d \n%d %d\n%d %d\n%d %d\n",batchCount, m,n,k, lda, tda, ldb,tdb, ldc,tdc);

  cublasDgemmBatched(handle,op_t1,op_t2,m,n,k,&alpha,(const double**) d_Aarray,lda, (const double**) d_Barray,ldb,&beta,(double**) d_Carray,ldc,batchCount);

  //printf("after dgemm\n");
  cudaDeviceSynchronize();
  
  //cudaFree(Aarray);
  //cudaFree(Barray);
  //cudaFree(Carray);
  
  //cudaFree(d_Aarray);
  //cudaFree(d_Barray);
  //cudaFree(d_Carray);
  //cublasDestroy(handle);
  
  
}

extern "C" void cublasDgemmBatched_finalize ()
{



  if (alreadyAllocated){
  
    cudaFree(Aarray);
    cudaFree(Barray);
    cudaFree(Carray);
    
    cudaFree(d_Aarray);
    cudaFree(d_Barray);
    cudaFree(d_Carray);
    cublasDestroy(handle);

  }
  
}
