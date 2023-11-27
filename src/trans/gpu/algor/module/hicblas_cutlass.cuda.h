#ifdef USE_CUTLASS
#include "hicblas.h"
#include "cutlass/gemm/device/gemm.h"

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
      OperatorClass, cutlass::arch::Sm50,                           //
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
  HIC_CHECK(cudaDeviceSynchronize());
}

void cutlass_sgemm_wrapper_grouped(char transa, char transb,
                                   int m, int *n, int *k, float alpha,
                                   const float *A, int lda, int *offsetsA,
                                   const float *B, int ldb, int *offsetsB, float beta,
                                   float *C, int ldc, int *offsetsC,
                                   int batchCount) {

  if (transa == 'N' && transb == 'N')
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_N, CUBLAS_OP_N>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
        batchCount);
  else if (transa == 'N' && transb == 'T')
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_N, CUBLAS_OP_T>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
        batchCount);
  else if (transa == 'T' && transb == 'N')
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_T, CUBLAS_OP_N>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
        batchCount);
  else if (transa == 'T' && transb == 'T')
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_T, CUBLAS_OP_T>(
        m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
        batchCount);
  else
    assert(false);
}

#endif
