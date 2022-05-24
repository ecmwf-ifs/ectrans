#include "cublas_v2.h"
#include <stdio.h>

#define CUDA_CHECK(e)                                                          \
  {                                                                            \
    cudaError_t err = (e);                                                     \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error: %s, line %d, %s: %s\n", __FILE__, __LINE__, \
              #e, cudaGetErrorString(err));                                    \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }
#define CUBLAS_CHECK(e)                                                        \
  {                                                                            \
    cublasStatus_t err = (e);                                                  \
    if (err != CUBLAS_STATUS_SUCCESS) {                                        \
      fprintf(stderr, "CUBLAS error: %s, line %d, %s: %i\n", __FILE__,         \
              __LINE__, #e, err);                                              \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

extern "C" {
void cublas_dgemm_wrapper(cublasOperation_t transa, cublasOperation_t transb,
                          int m, int n, int k, double alpha, const double *A,
                          int lda, int tda, const double *B, int ldb, int tdb,
                          double beta, double *C, int ldc, int tdc,
                          int batchCount) {
  static cublasHandle_t handle = nullptr;
  if (!handle)
    CUBLAS_CHECK(cublasCreate(&handle));

  CUBLAS_CHECK(cublasDgemmStridedBatched(handle, transa, transb, m, n, k,
                                         &alpha, A, lda, tda, B, ldb, tdb,
                                         &beta, C, ldc, tdc, batchCount));
  CUDA_CHECK(cudaDeviceSynchronize());
}

void cublas_sgemm_wrapper(cublasOperation_t transa, cublasOperation_t transb,
                          int m, int n, int k, float alpha, const float *A,
                          int lda, int tda, const float *B, int ldb, int tdb,
                          float beta, float *C, int ldc, int tdc,
                          int batchCount) {
  static cublasHandle_t handle = nullptr;
  if (!handle)
    CUBLAS_CHECK(cublasCreate(&handle));

  CUBLAS_CHECK(cublasSgemmStridedBatched(handle, transa, transb, m, n, k,
                                         &alpha, A, lda, tda, B, ldb, tdb,
                                         &beta, C, ldc, tdc, batchCount));
  CUDA_CHECK(cudaDeviceSynchronize());
}
}
