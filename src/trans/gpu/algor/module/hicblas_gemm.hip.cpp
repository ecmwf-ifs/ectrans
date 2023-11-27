// (C) Copyright 2000- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdio.h>
#include <stdlib.h>

#include "hicblas.h"

#ifdef USE_CUTLASS
#include "hicblas_cutlass.cuda.h"
constexpr bool use_cutlass = false;
#else
constexpr bool use_cutlass = false;
#endif

bool hip_alreadyAllocated_sgemm=false;
bool hip_alreadyAllocated_sgemm_handle=false;

bool hip_alreadyAllocated_dsgemm=false;
bool hip_alreadyAllocated_dgemm_handle=false;

hipblasHandle_t handle_hip_sgemm;
hipblasHandle_t handle_hip_dgemm;


extern "C" void hipblas_sgemm_wrapper (char transa, char transb,
                                       int m, int n,int k, float alpha,
                                       const float *A, int lda, int tda,
                                       const float *B, int ldb, int tdb, float beta,
                                       float *C, int ldc, int tdc, int batchCount)
{


  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=HIPBLAS_OP_T;

  if (!hip_alreadyAllocated_sgemm_handle){
    hipblasCreate(&handle_hip_sgemm);
    hip_alreadyAllocated_sgemm_handle=true;
  }
  hipblasSgemmStridedBatched(handle_hip_sgemm,op_t1,op_t2,m,n,k,
                             &alpha,(const float *) A,lda,tda, (const float *) B,ldb,tdb,
                             &beta,(float*) C,ldc,tdc,batchCount);

}

extern "C" void hipblas_dgemm_wrapper (char transa, char transb,
                                       int m, int n,int k, double alpha,
                                       const double *A, int lda, int tda,
                                       const double *B, int ldb, int tdb, double beta,
                                       double *C, int ldc, int tdc, int batchCount)
{


  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=HIPBLAS_OP_T;

  if (!hip_alreadyAllocated_dgemm_handle){
    hipblasCreate(&handle_hip_dgemm);
    hip_alreadyAllocated_dgemm_handle=true;
  }
  hipblasDgemmStridedBatched(handle_hip_sgemm,op_t1,op_t2,m,n,k,
                             &alpha,(const double *) A,lda,tda, (const double *) B,ldb,tdb,
                             &beta,(double *) C,ldc,tdc,batchCount);

}


extern "C" void hipblasDgemmGrouped_wrapper(char transa, char transb,
                                  int m, int *n, int *k,
                                  double alpha,
                                  const double *A, int lda, int *offsetsA,
                                  const double *B, int ldb, int *offsetsB,
                                  double beta,
                                  double *C, int ldc, int *offsetsC,
                                  int batchCount) {

  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=HIPBLAS_OP_T;

  static hipblasHandle_t handle_dgemm_grouped = nullptr;
  if (!handle_dgemm_grouped)
    HICBLAS_CHECK(hipblasCreate(&handle_dgemm_grouped));

  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    HICBLAS_CHECK(hipblasDgemm(handle_dgemm_grouped, op_t1, op_t2, m, n[i], k[i], &alpha,
                             A + offsetsA[i], lda, B + offsetsB[i], ldb, &beta,
                             C + offsetsC[i], ldc));
  }
}


extern "C" void hipblasSgemmGrouped_wrapper(char transa, char transb,
                                  int m, int *n, int *k, float alpha,
                                  const float *A, int lda, int *offsetsA,
                                  const float *B, int ldb, int *offsetsB, float beta,
                                  float *C, int ldc, int *offsetsC,
                                  int batchCount) {

  hipblasOperation_t op_t1=HIPBLAS_OP_N, op_t2=HIPBLAS_OP_N;

  if (transa=='T' || transa=='t')
    op_t1=HIPBLAS_OP_T;

  if (transb=='T' || transb=='t')
    op_t2=HIPBLAS_OP_T;

  static hipblasHandle_t handle_sgemm_grouped = nullptr;
  if (!handle_sgemm_grouped)
    HICBLAS_CHECK(hipblasCreate(&handle_sgemm_grouped));

  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    HICBLAS_CHECK(hipblasSgemm(handle_sgemm_grouped, op_t1, op_t2, m, n[i], k[i], &alpha,
                             A + offsetsA[i], lda, B + offsetsB[i], ldb, &beta,
                             C + offsetsC[i], ldc));
  }
}

extern "C" void blas_dgemm_wrapper_grouped(char transa, char transb,
                                int m, int *n, int *k, double alpha,
                                const double *A, int lda, int *offsetsA,
                                const double *B, int ldb, int *offsetsB, double beta,
                                double *C, int ldc, int *offsetsC,
                                int batchCount) {
    hipblasDgemmGrouped_wrapper(transa, transb, m, n, k, alpha, A, lda, offsetsA, B,
                                ldb, offsetsB, beta, C, ldc, offsetsC, batchCount);
}

extern "C" void blas_sgemm_wrapper_grouped(char transa, char transb,
                                int m, int *n, int *k, float alpha,
                                const float *A, int lda, int *offsetsA,
                                const float *B, int ldb, int *offsetsB, float beta,
                                float *C, int ldc, int *offsetsC,
                                int batchCount) {
#ifdef USE_CUTLASS
    cutlass_sgemm_wrapper_grouped(transa, transb, m, n, k, alpha, A, lda, offsetsA,
                                  B, ldb, offsetsB, beta, C, ldc, offsetsC, batchCount);
#else
    hipblasSgemmGrouped_wrapper(transa, transb, m, n, k, alpha, A, lda, offsetsA, B,
                                ldb, offsetsB, beta, C, ldc, offsetsC, batchCount);
#endif
}



extern "C" void hipblasSgemmBatched_finalize ()
{

#ifdef FALSE
  if (hip_alreadyAllocated_sgemm){

    hipFree(Aarray_sgemm_hip);
    hipFree(Barray_sgemm_hip);
    hipFree(Carray_sgemm_hip);

    hipFree(d_Aarray_sgemm_hip);
    hipFree(d_Barray_sgemm_hip);
    hipFree(d_Carray_sgemm_hip);

  }
#endif

  if (hip_alreadyAllocated_sgemm_handle){
    hipblasDestroy(handle_hip_sgemm);
  }

}


