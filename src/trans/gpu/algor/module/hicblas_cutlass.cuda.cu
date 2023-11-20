#ifdef USE_CUTLASS
#include "hicblas.h"
#include "cutlass/gemm/device/gemm.h"



template <cublasOperation_t TransA, cublasOperation_t TransB>
void cutlass_sgemm_wrapper_grouped_v(int m, int *n, int *k, float alpha,
                                     const float *A, int lda, int tda,
                                     const float *B, int ldb, int tdb,
                                     float beta, float *C, int ldc, int tdc,
                                     int batchCount) {
  // this was verified using Volta and uses FP32
  constexpr int AlignmentA = 1;
  constexpr int AlignmentB = 1;x
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
    CUTLASS_CHECK(gemm_op({//
                           {(m + sz_align - 1) / sz_align * sz_align,
                            (n[i] + sz_align - 1) / sz_align * sz_align,
                            (k[i] + sz_align - 1) / sz_align * sz_align},
                           {const_cast<float *>(A + i * tda), lda},
                           {const_cast<float *>(B + i * tdb), ldb},
                           {C + i * tdc, ldc},
                           {C + i * tdc, ldc},
                           {alpha, beta}}));
  }
  CUDA_CHECK(cudaDeviceSynchronize());
}
void cutlass_sgemm_wrapper_grouped(cublasOperation_t transa,
                                   cublasOperation_t transb, int m, int *n,
                                   int *k, float alpha, const float *A, int lda,
                                   int tda, const float *B, int ldb, int tdb,
                                   float beta, float *C, int ldc, int tdc,
                                   int batchCount) {
  if (transa == CUBLAS_OP_N && transb == CUBLAS_OP_N)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_N, CUBLAS_OP_N>(
        m, n, k, alpha, A, lda, tda, B, ldb, tdb, beta, C, ldc, tdc,
        batchCount);
  else if (transa == CUBLAS_OP_N && transb == CUBLAS_OP_T)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_N, CUBLAS_OP_T>(
        m, n, k, alpha, A, lda, tda, B, ldb, tdb, beta, C, ldc, tdc,
        batchCount);
  else if (transa == CUBLAS_OP_T && transb == CUBLAS_OP_N)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_T, CUBLAS_OP_N>(
        m, n, k, alpha, A, lda, tda, B, ldb, tdb, beta, C, ldc, tdc,
        batchCount);
  else if (transa == CUBLAS_OP_T && transb == CUBLAS_OP_T)
    cutlass_sgemm_wrapper_grouped_v<CUBLAS_OP_T, CUBLAS_OP_T>(
        m, n, k, alpha, A, lda, tda, B, ldb, tdb, beta, C, ldc, tdc,
        batchCount);
  else
    assert(false);
}

#endif