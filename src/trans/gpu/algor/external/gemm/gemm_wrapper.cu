#include "cublas_v2.h"
#include "cutlass/gemm/device/gemm.h"
#include <stdio.h>

constexpr bool use_cutlass = true;

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
#define CUTLASS_CHECK(e)                                                       \
  {                                                                            \
    cutlass::Status err = (e);                                                 \
    if (err != cutlass::Status::kSuccess) {                                    \
      fprintf(stderr, "CUTLASS error: %s, line %d, %s: %i\n", __FILE__,        \
              __LINE__, #e, (int)err);                                         \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

template <cublasOperation_t TransA, cublasOperation_t TransB>
void cutlass_sgemm_wrapper_grouped_v(int m, int *n, int *k, float alpha,
                                     const float *A, int lda, int *offsetsA,
                                     const float *B, int ldb, int *offsetsB,
                                     float beta, float *C, int ldc,
                                     int *offsetsC, int batchCount) {
#if 0
  // we will enable this later (this ifdefs did not work, so I am going to enable this properly ltaer)
  // this was verified using Ampere and uses 3XTF32
  constexpr int AlignmentA = 4;
  constexpr int AlignmentB = 4;
  using ThreadblockShape = cutlass::gemm::GemmShape<128, 64, 32>;
  using WarpShape = cutlass::gemm::GemmShape<64, 32, 32>;
  using InstructionShape = cutlass::gemm::GemmShape<16, 8, 8>;
  using OperatorClass = cutlass::arch::OpClassTensorOp;
  using MyOp = cutlass::arch::OpMultiplyAddFastF32;

  using Gemm = cutlass::gemm::device::Gemm<
      float,
      std::conditional_t<TransA == CUBLAS_OP_N, cutlass::layout::ColumnMajor,
                         cutlass::layout::RowMajor>, //
      float,
      std::conditional_t<TransB == CUBLAS_OP_N, cutlass::layout::ColumnMajor,
                         cutlass::layout::RowMajor>, //
      float, cutlass::layout::ColumnMajor,           //
      float,                                         //
      OperatorClass, cutlass::arch::Sm80,            //
      ThreadblockShape, WarpShape, InstructionShape, //
      cutlass::epilogue::thread::LinearCombination<  //
          float,                                     //
          128 / cutlass::sizeof_bits<float>::value,
          float,                                                    //
          float                                                     //
          >,                                                        //
      cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>, //
      3,                                                            //
      AlignmentA,                                                   //
      AlignmentB,                                                   //
      true,                                                         //
      MyOp                                                          //
      >;
  constexpr int sz_align = 8;
#else
  // this was verified using Volta and uses FP32
  constexpr int AlignmentA = 1;
  constexpr int AlignmentB = 1;
  using ThreadblockShape = cutlass::gemm::GemmShape<128, 128, 8>;
  using WarpShape = cutlass::gemm::GemmShape<32, 32, 8>;
  using InstructionShape = cutlass::gemm::GemmShape<1, 1, 1>;
  using OperatorClass = cutlass::arch::OpClassSimt;
  using MyOp = cutlass::arch::OpMultiplyAdd;

  using Gemm = cutlass::gemm::device::Gemm<
      float, //
      std::conditional_t<TransA == CUBLAS_OP_N, cutlass::layout::ColumnMajor,
                         cutlass::layout::RowMajor>, //
      float,                                         //
      std::conditional_t<TransB == CUBLAS_OP_N, cutlass::layout::ColumnMajor,
                         cutlass::layout::RowMajor>,                //
      float, cutlass::layout::ColumnMajor,                          //
      float,                                                        //
      OperatorClass, cutlass::arch::Sm70,                           //
      ThreadblockShape, WarpShape, InstructionShape,                //
      cutlass::epilogue::thread::LinearCombination<                 //
          float,                                                    //
          1,                                                        //
          float,                                                    //
          float                                                     //
          >,                                                        //
      cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>, //
      2,                                                            //
      AlignmentA,                                                   //
      AlignmentB,                                                   //
      true,                                                         //
      MyOp                                                          //
      >;
  constexpr int sz_align = 1;
#endif

  Gemm gemm_op;
  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    CUTLASS_CHECK(gemm_op({//
                           {(m + sz_align - 1) / sz_align * sz_align,
                            (n[i] + sz_align - 1) / sz_align * sz_align,
                            (k[i] + sz_align - 1) / sz_align * sz_align},
                           {const_cast<float *>(A + offsetsA[i]), lda},
                           {const_cast<float *>(B + offsetsB[i]), ldb},
                           {C + offsetsC[i], ldc},
                           {C + offsetsC[i], ldc},
                           {alpha, beta}}));
  }
  CUDA_CHECK(cudaDeviceSynchronize());
}
void cutlass_sgemm_wrapper_grouped(cublasOperation_t transa,
                                   cublasOperation_t transb, int m, int *n,
                                   int *k, float alpha, const float *A, int lda,
                                   int *offsetsA, const float *B, int ldb,
                                   int *offsetsB, float beta, float *C, int ldc,
                                   int *offsetsC, int batchCount) {
  if (transa == CUBLAS_OP_N && transb == CUBLAS_OP_N)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_N, CUBLAS_OP_N>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc,
        offsetsC, batchCount);
  else if (transa == CUBLAS_OP_N && transb == CUBLAS_OP_T)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_N, CUBLAS_OP_T>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc,
        offsetsC, batchCount);
  else if (transa == CUBLAS_OP_T && transb == CUBLAS_OP_N)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_T, CUBLAS_OP_N>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc,
        offsetsC, batchCount);
  else if (transa == CUBLAS_OP_T && transb == CUBLAS_OP_T)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_T, CUBLAS_OP_T>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc,
        offsetsC, batchCount);
  else
    assert(false);
}
void cublas_sgemm_wrapper_grouped(cublasOperation_t transa,
                                  cublasOperation_t transb, int m, int *n,
                                  int *k, float alpha, const float *A, int lda,
                                  int *offsetsA, const float *B, int ldb,
                                  int *offsetsB, float beta, float *C, int ldc,
                                  int *offsetsC, int batchCount) {
  static cublasHandle_t handle = nullptr;
  if (!handle)
    CUBLAS_CHECK(cublasCreate(&handle));

  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    CUBLAS_CHECK(cublasSgemm(handle, transa, transb, m, n[i], k[i], &alpha,
                             A + offsetsA[i], lda, B + offsetsB[i], ldb, &beta,
                             C + offsetsC[i], ldc));
  }
}
void cublas_dgemm_wrapper_grouped(cublasOperation_t transa,
                                  cublasOperation_t transb, int m, int *n,
                                  int *k, double alpha, const double *A,
                                  int lda, int *offsetsA, const double *B,
                                  int ldb, int *offsetsB, double beta,
                                  double *C, int ldc, int *offsetsC,
                                  int batchCount) {
  static cublasHandle_t handle = nullptr;
  if (!handle)
    CUBLAS_CHECK(cublasCreate(&handle));

  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    CUBLAS_CHECK(cublasDgemm(handle, transa, transb, m, n[i], k[i], &alpha,
                             A + offsetsA[i], lda, B + offsetsB[i], ldb, &beta,
                             C + offsetsC[i], ldc));
  }
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
}

void blas_sgemm_wrapper_grouped(cublasOperation_t transa,
                                cublasOperation_t transb, int m, int *n, int *k,
                                float alpha, const float *A, int lda,
                                int *offsetsA, const float *B, int ldb,
                                int *offsetsB, float beta, float *C, int ldc,
                                int *offsetsC, int batchCount) {
  if (use_cutlass)
    cutlass_sgemm_wrapper_grouped(transa, transb, m, n, k, alpha, A, lda,
                                  offsetsA, B, ldb, offsetsB, beta, C, ldc,
                                  offsetsC, batchCount);
  else
    cublas_sgemm_wrapper_grouped(transa, transb, m, n, k, alpha, A, lda,
                                 offsetsA, B, ldb, offsetsB, beta, C, ldc,
                                 offsetsC, batchCount);
}
void blas_dgemm_wrapper_grouped(cublasOperation_t transa,
                                cublasOperation_t transb, int m, int *n, int *k,
                                double alpha, const double *A, int lda,
                                int *offsetsA, const double *B, int ldb,
                                int *offsetsB, double beta, double *C, int ldc,
                                int *offsetsC, int batchCount) {
  cublas_dgemm_wrapper_grouped(transa, transb, m, n, k, alpha, A, lda, offsetsA,
                               B, ldb, offsetsB, beta, C, ldc, offsetsC,
                               batchCount);
}
}
