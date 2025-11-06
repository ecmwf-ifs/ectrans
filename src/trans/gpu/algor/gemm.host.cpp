// (C) Copyright 2000- ECMWF.
// (C) Copyright 2024- NVIDIA.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.


#include <stdlib.h>

extern "C" void dgemm_(
  const char *transa, const char *transb,
  const int *m, const int *n, const int *k,
  double *alpha,
  const double *A, const int *lda,
  const double *B, const int *ldb,
  double *beta,
  double *C, const int *ldc
);

extern "C" void sgemm_(
  const char *transa, const char *transb,
  const int *m, const int *n, const int *k,
  float *alpha,
  const float *A, const int *lda,
  const float *B, const int *ldb,
  float *beta,
  float *C, const int *ldc
);

extern "C" {
void hipblas_dgemm_wrapper(char transa, char transb, int m, int n, int k,
                           double alpha, const double *A, int lda, int tda,
                           const double *B, int ldb, int tdb, double beta,
                           double *C, int ldc, int tdc, int batchCount,
                           size_t stream, void *growing_allocator) {
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

void hipblas_sgemm_wrapper(char transa, char transb, int m, int n, int k,
                           float alpha, const float *A, int lda, int tda,
                           const float *B, int ldb, int tdb, float beta,
                           float *C, int ldc, int tdc, int batchCount,
                           void *growing_allocator) {
  sgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

void hipblas_sgemm_wrapper_grouped(
    int resol_id, int blas_id, char transa, char transb, int m, const int *n,
    const int *k, float alpha, const float *A, int lda, const int64_t *offsetsA,
    const float *B, const int *ldb, const int64_t *offsetsB, float beta,
    float *C, int ldc, const int64_t *offsetsC, int batchCount, size_t stream,
    void *growing_allocator) {

  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    sgemm_(&transa, &transb, &m, &n[i], &k[i], &alpha, A + offsetsA[i], &lda, B + offsetsB[i],
         &ldb[i], &beta, C + offsetsC[i], &ldc);
  }
}

void hipblas_dgemm_wrapper_grouped(int resol_id, int blas_id, char transa,
                                   char transb, int m, const int *n,
                                   const int *k, double alpha, const double *A,
                                   int lda, const int64_t *offsetsA,
                                   const double *B, const int *ldb,
                                   const int64_t *offsetsB, double beta,
                                   double *C, int ldc, const int64_t *offsetsC,
                                   int batchCount, size_t stream,
                                   void *growing_allocator) {
  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    dgemm_(&transa, &transb, &m, &n[i], &k[i], &alpha, A + offsetsA[i], &lda, B + offsetsB[i],
         &ldb[i], &beta, C + offsetsC[i], &ldc);
  }
}

void clean_gemm(int resol_id) {
  // Just here to satisfy compilation -> there's nothing to clean on the host version
}
}
