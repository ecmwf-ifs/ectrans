#ifdef USE_CUTLASS
//#include "hicblas.h"
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

#ifdef USE_CUTLASS_3XTF32
constexpr bool use_3xtf32 = true;
#else
constexpr bool use_3xtf32 = false;
#endif


template <typename CutlassGemm>
CutlassGemm &get_cutlass_handle() {
  static auto handle = std::make_unique<CutlassGemm>();
  return *handle;
}

namespace detail {

enum class CutlassType { cutlass_3xtf32, cutlass_fp32 };

template <CutlassType, cublasOperation_t TransA, cublasOperation_t TransB>
class cutlass_sgemm_grouped;

template <cublasOperation_t TransA, cublasOperation_t TransB>
class cutlass_sgemm_grouped<CutlassType::cutlass_3xtf32, TransA, TransB> {
  // this was verified using Ampere and uses 3XTF32
  static constexpr int AlignmentA = 4;
  static constexpr int AlignmentB = 4;
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
  static constexpr int sz_align = 8;

public:
  void operator()(cudaStream_t stream, int m, int n, int k, float alpha,
                  const float *A, int lda, const float *B, int ldb, float beta,
                  float *C, int ldc) const {
    auto &gemm_op = get_cutlass_handle<Gemm>();
    CUTLASS_CHECK(gemm_op(
        {//
         {(m + sz_align - 1) / sz_align * sz_align,
          (n + sz_align - 1) / sz_align * sz_align,
          (k + sz_align - 1) / sz_align * sz_align},
         {const_cast<float *>(A), lda},
         {const_cast<float *>(B), ldb},
         {C, ldc},
         {C, ldc},
         {alpha, beta}},
        nullptr, stream));
  }
};
template <cublasOperation_t TransA, cublasOperation_t TransB>
class cutlass_sgemm_grouped<CutlassType::cutlass_fp32, TransA, TransB> {
  // this was verified using Volta and uses FP32
  static constexpr int AlignmentA = 1;
  static constexpr int AlignmentB = 1;
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
  static constexpr int sz_align = 1;

public:
  void operator()(cudaStream_t stream, int m, int n, int k, float alpha,
                  const float *A, int lda, const float *B, int ldb, float beta,
                  float *C, int ldc) const {
    auto &gemm_op = get_cutlass_handle<Gemm>();
    CUTLASS_CHECK(gemm_op(
        {//
         {(m + sz_align - 1) / sz_align * sz_align,
          (n + sz_align - 1) / sz_align * sz_align,
          (k + sz_align - 1) / sz_align * sz_align},
         {const_cast<float *>(A), lda},
         {const_cast<float *>(B), ldb},
         {C, ldc},
         {C, ldc},
         {alpha, beta}},
        nullptr, stream));
  }
};

} // namespace detail
template <cublasOperation_t TransA, cublasOperation_t TransB>
void cutlass_sgemm_wrapper_grouped_op(int blas_id, int m, int *n, int *k,
                                      float alpha, const float *A, int lda,
                                      int *offsetsA, const float *B, int ldb,
                                      int *offsetsB, float beta, float *C,
                                      int ldc, int *offsetsC, int batchCount,
                                      cudaStream_t stream,
                                      void *growing_allocator) {
  using namespace detail;
  int device;
  HIC_CHECK(cudaGetDevice(&device));
  int capability_major;
  HIC_CHECK(cudaDeviceGetAttribute(&capability_major,
                                    cudaDevAttrComputeCapabilityMajor, device));
  if (capability_major >= 8 && use_3xtf32)
    run_group_graph(cutlass_sgemm_grouped<detail::CutlassType::cutlass_3xtf32,
                                          TransA, TransB>(),
                    m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
                    ldc, offsetsC, batchCount, stream, blas_id,
                    growing_allocator);
  else
    run_group_graph(cutlass_sgemm_grouped<detail::CutlassType::cutlass_fp32,
                                          TransA, TransB>(),
                    m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
                    ldc, offsetsC, batchCount, stream, blas_id,
                    growing_allocator);
}

void cutlass_sgemm_wrapper_grouped(int blas_id, char transa, char transb,
                                   int m, int *n, int *k, float alpha,
                                   const float *A, int lda, int *offsetsA,
                                   const float *B, int ldb, int *offsetsB, float beta,
                                   float *C, int ldc, int *offsetsC,
                                   int batchCount, cudaStream_t stream,
                                   void *growing_allocator) {

  if (transa == 'N' && transb == 'N')
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_N, CUBLAS_OP_N>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount, stream, growing_allocator);
  else if (transa == 'N' && transb == 'T')
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_N, CUBLAS_OP_T>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount, stream, growing_allocator);
  else if (transa == 'T' && transb == 'N')
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_T, CUBLAS_OP_N>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount, stream, growing_allocator);
  else if (transa == 'T' && transb == 'T')
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_T, CUBLAS_OP_T>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount, stream, growing_allocator);
  else
    assert(false);
}
//}

#endif
