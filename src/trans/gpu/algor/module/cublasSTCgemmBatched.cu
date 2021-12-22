//
// Wrapper for cublasSTCgemm function. 
//
// Sam Hatfield, ECMWF
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "cublas_v2.h" 

bool alreadyAllocated_stcgemm = false;
bool alreadyAllocated_stcgemm_handle = false;

half **d_Aarray_stcgemm;
half **d_Barray_stcgemm;
float **d_Carray_stcgemm;

half **Aarray_stcgemm;
half **Barray_stcgemm;
float **Carray_stcgemm;

cublasHandle_t handle_stcgemm;

extern "C" void cublasSTCgemmBatched_wrapper(
        char transa, char transb,
        int m, int n, int k,
        float alpha,
        const half *A, int lda, int tda,
        const half *B, int ldb, int tdb,
        float beta,
        float *C, int ldc, int tdc,
        int batchCount
){
  // Define CUBLAS operation handles
  cublasOperation_t op_t1, op_t2;

  // Decide whether to transpose matrices or not
  op_t1 = (transa == 'T' || transa == 't') ? CUBLAS_OP_T : CUBLAS_OP_N;
  op_t2 = (transb == 'T' || transb == 't') ? CUBLAS_OP_T : CUBLAS_OP_N;

  // Initialize CUBLAS handle
  if (!alreadyAllocated_stcgemm_handle) {
    cublasCreate(&handle_stcgemm);
    alreadyAllocated_stcgemm_handle = true;
  }

  // Allocate host arrays
  if (!alreadyAllocated_stcgemm) {
    cudaMallocHost(&Aarray_stcgemm,batchCount*sizeof(half*));
    cudaMallocHost(&Barray_stcgemm,batchCount*sizeof(half*));
    cudaMallocHost(&Carray_stcgemm,batchCount*sizeof(float*));
    alreadyAllocated_stcgemm = true;
  }

  // Allocate device arrays
  cudaMalloc(&d_Aarray_stcgemm, batchCount*sizeof(half*));
  cudaMalloc(&d_Barray_stcgemm, batchCount*sizeof(half*));
  cudaMalloc(&d_Carray_stcgemm, batchCount*sizeof(float*));

  // Transfer data from input arrays to host arrays
  for (int i = 0; i < batchCount; i++) {
    Aarray_stcgemm[i] = (half*) &(A[i*lda*tda]);
    Barray_stcgemm[i] = (half*) &(B[i*ldb*tdb]);
    Carray_stcgemm[i] = (float*) &(C[i*ldc*tdc]);
  }

  // Transfer data from host arrays to device arrays
  cudaMemcpy(d_Aarray_stcgemm, Aarray_stcgemm, batchCount*sizeof(half*), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Barray_stcgemm, Barray_stcgemm, batchCount*sizeof(half*), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Carray_stcgemm, Carray_stcgemm, batchCount*sizeof(float*), cudaMemcpyHostToDevice);

  // Perform batched SGEMM
  cublasGemmBatchedEx(handle_stcgemm,
    op_t1, op_t2,
    m, n, k,
    (const void*)&alpha,
    (const void**)d_Aarray_stcgemm, CUDA_R_16F, lda,
    (const void**)d_Barray_stcgemm, CUDA_R_16F, ldb,
    (const void*)&beta,
    (void**)d_Carray_stcgemm, CUDA_R_32F, ldc,
    batchCount,
    CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT_TENSOR_OP
  );

  cudaDeviceSynchronize();
  
  // Free device arrays
  cudaFree(d_Aarray_stcgemm);
  cudaFree(d_Barray_stcgemm);
  cudaFree(d_Carray_stcgemm);
}

extern "C" void cublasSTCgemmBatched_finalize() {
  if (alreadyAllocated_stcgemm) {
    cudaFree(Aarray_stcgemm);
    cudaFree(Barray_stcgemm);
    cudaFree(Carray_stcgemm);
    
    cudaFree(d_Aarray_stcgemm);
    cudaFree(d_Barray_stcgemm);
    cudaFree(d_Carray_stcgemm);
  }

  if (alreadyAllocated_stcgemm_handle) {
    cublasDestroy(handle_stcgemm);
  }
}
