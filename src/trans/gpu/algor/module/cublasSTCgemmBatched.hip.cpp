//
// Wrapper for cublasSTCgemm function. 
//
// Sam Hatfield, ECMWF
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "hipblas.h" 

bool alreadyAllocated_stcgemm = false;
bool alreadyAllocated_stcgemm_handle = false;

half **d_Aarray_stcgemm;
half **d_Barray_stcgemm;
float **d_Carray_stcgemm;

half **Aarray_stcgemm;
half **Barray_stcgemm;
float **Carray_stcgemm;

hipblasHandle_t handle_stcgemm;

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
  hipblasOperation_t op_t1, op_t2;

  // Decide whether to transpose matrices or not
  op_t1 = (transa == 'T' || transa == 't') ? HIPBLAS_OP_T : HIPBLAS_OP_N;
  op_t2 = (transb == 'T' || transb == 't') ? HIPBLAS_OP_T : HIPBLAS_OP_N;

  // Initialize CUBLAS handle
  if (!alreadyAllocated_stcgemm_handle) {
    hipblasCreate(&handle_stcgemm);
    alreadyAllocated_stcgemm_handle = true;
  }

  // Allocate host arrays
  if (!alreadyAllocated_stcgemm) {
    hipHostMalloc(&Aarray_stcgemm,batchCount*sizeof(half*));
    hipHostMalloc(&Barray_stcgemm,batchCount*sizeof(half*));
    hipHostMalloc(&Carray_stcgemm,batchCount*sizeof(float*));
    alreadyAllocated_stcgemm = true;
  }

  // Allocate device arrays
  hipMalloc(&d_Aarray_stcgemm, batchCount*sizeof(half*));
  hipMalloc(&d_Barray_stcgemm, batchCount*sizeof(half*));
  hipMalloc(&d_Carray_stcgemm, batchCount*sizeof(float*));

  // Transfer data from input arrays to host arrays
  for (int i = 0; i < batchCount; i++) {
    Aarray_stcgemm[i] = (half*) &(A[i*lda*tda]);
    Barray_stcgemm[i] = (half*) &(B[i*ldb*tdb]);
    Carray_stcgemm[i] = (float*) &(C[i*ldc*tdc]);
  }

  // Transfer data from host arrays to device arrays
  hipMemcpy(d_Aarray_stcgemm, Aarray_stcgemm, batchCount*sizeof(half*), hipMemcpyHostToDevice);
  hipMemcpy(d_Barray_stcgemm, Barray_stcgemm, batchCount*sizeof(half*), hipMemcpyHostToDevice);
  hipMemcpy(d_Carray_stcgemm, Carray_stcgemm, batchCount*sizeof(float*), hipMemcpyHostToDevice);

  // Perform batched SGEMM
  hipblasGemmBatchedEx(handle_stcgemm,
    op_t1, op_t2,
    m, n, k,
    (const void*)&alpha,
    (const void**)d_Aarray_stcgemm, HIPBLAS_R_16F, lda,
    (const void**)d_Barray_stcgemm, HIPBLAS_R_16F, ldb,
    (const void*)&beta,
    (void**)d_Carray_stcgemm, HIPBLAS_R_32F, ldc,
    batchCount,
    CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT_TENSOR_OP
  );

  hipDeviceSynchronize();
  
  // Free device arrays
  hipFree(d_Aarray_stcgemm);
  hipFree(d_Barray_stcgemm);
  hipFree(d_Carray_stcgemm);
}

extern "C" void cublasSTCgemmBatched_finalize() {
  if (alreadyAllocated_stcgemm) {
    hipFree(Aarray_stcgemm);
    hipFree(Barray_stcgemm);
    hipFree(Carray_stcgemm);
    
    hipFree(d_Aarray_stcgemm);
    hipFree(d_Barray_stcgemm);
    hipFree(d_Carray_stcgemm);
  }

  if (alreadyAllocated_stcgemm_handle) {
    hipblasDestroy(handle_stcgemm);
  }
}
