//
// Wrapper for Tensor Core GEMM function.
//
// Sam Hatfield, ECMWF
// Alan Gray, NVIDIA
//

#include <stdio.h>
#include "cublas_v2.h" 


bool alreadyAllocated_tcgemm=false;
bool alreadyAllocated_tcgemm_handle=false;

// Device arrays
//half **d_Aarray_h;
//half **d_Barray_h;
float **d_Aarray_h;
float **d_Barray_h;
float **d_Aarray_tcgemm;
float **d_Barray_tcgemm;
float **d_Carray_tcgemm;

// Host arrays
float **Aarray_tcgemm;
float **Barray_tcgemm;
float **Carray_tcgemm;

cublasHandle_t handle_tcgemm;

// Converts from single-precision to half-precision (CUDA kernel)
__global__ void float2half(half *out, const float *in, int n) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < n) {
    out[idx] = __float2half(in[idx]);
  }
}



extern "C" void cublasTCgemmBatched_wrapper(char transa, char transb,
                                            int m, int n, int k,
                                            float alpha,
                                            const float *A, int lda, int tda,
                                            const float *B, int ldb, int tdb,
                                            float beta,
                                            float *C, int ldc, int tdc,
                                            int batchCount)
{
  fprintf(stderr, "Using Tensor Core\n");

  // Set transpose operation parameters
  cublasOperation_t op_t1 = (transa == 'T' || transa == 't') ? CUBLAS_OP_T : CUBLAS_OP_N;
  cublasOperation_t op_t2 = (transb == 'T' || transb == 't') ? CUBLAS_OP_T : CUBLAS_OP_N;

  if (!alreadyAllocated_tcgemm_handle) {
    cublasCreate(&handle_tcgemm);
    alreadyAllocated_tcgemm_handle=true;
  }

  //cublasSetMathMode(handle_tcgemm, CUBLAS_TENSOR_OP_MATH);

  if (!alreadyAllocated_tcgemm) {
    // Allocate host arrays specifically for host->device transfer
    cudaMallocHost(&Aarray_tcgemm, batchCount*sizeof(float*));
    cudaMallocHost(&Barray_tcgemm, batchCount*sizeof(float*));
    cudaMallocHost(&Carray_tcgemm, batchCount*sizeof(float*));
    alreadyAllocated_tcgemm=true;
  }

  // Allocate device arrays
  cudaMalloc(&d_Aarray_h, batchCount*sizeof(float*));
  cudaMalloc(&d_Barray_h, batchCount*sizeof(float*));
  cudaMalloc(&d_Aarray_tcgemm, batchCount*sizeof(float*));
  cudaMalloc(&d_Barray_tcgemm, batchCount*sizeof(float*));
  cudaMalloc(&d_Carray_tcgemm, batchCount*sizeof(float*));
 
  // Copy data from dummy arrays to host arrays
  for (int i = 0; i < batchCount; i++) {
    Aarray_tcgemm[i] = (float*) &(A[i*lda*tda]);
    Barray_tcgemm[i] = (float*) &(B[i*ldb*tdb]);
    Carray_tcgemm[i] = (float*) &(C[i*ldc*tdc]);
  }

  // Transfer arrays from host to device
  cudaMemcpy(d_Aarray_tcgemm, Aarray_tcgemm, batchCount*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Barray_tcgemm, Barray_tcgemm, batchCount*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Carray_tcgemm, Carray_tcgemm, batchCount*sizeof(float*), cudaMemcpyHostToDevice);

//  // Convert arrays to half-precision
//  for (int i = 0; i < batchCount; i++) {
//    float2half<<<(int)(m*k/256) + 1, 256 >>>(d_Aarray_h[i], d_Aarray_tcgemm[i], batchCount);
//    float2half<<<(int)(k*n/256) + 1, 256 >>>(d_Barray_h[i], d_Barray_tcgemm[i], batchCount);
//  }

  // Perform Tensor Core batched GEMM
  cublasGemmBatchedEx(handle_tcgemm, op_t1, op_t2,
                      m, n, k,
                      (const void *)&alpha,
                      (const void **)d_Aarray_h, CUDA_R_32F, lda,
                      (const void **)d_Barray_h, CUDA_R_32F, ldb,
                      (const void *)&beta,
                      (void **)d_Carray_tcgemm, CUDA_R_32F, ldc,
                      batchCount,
                      CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT);

  cudaDeviceSynchronize();

  cudaFree(d_Aarray_h);
  cudaFree(d_Barray_h);
  cudaFree(d_Aarray_tcgemm);
  cudaFree(d_Barray_tcgemm);
  cudaFree(d_Carray_tcgemm);
}

extern "C" void cublasTCgemmBatched_finalize()
{
  if (alreadyAllocated_tcgemm) {
    cudaFree(Aarray_tcgemm);
    cudaFree(Barray_tcgemm);
    cudaFree(Carray_tcgemm);
    
    cudaFree(d_Aarray_h);
    cudaFree(d_Barray_h);
    cudaFree(d_Aarray_tcgemm);
    cudaFree(d_Barray_tcgemm);
    cudaFree(d_Carray_tcgemm);
  }

  if (alreadyAllocated_tcgemm_handle) {
      cublasDestroy(handle_tcgemm);
  }
}

