#include <stdio.h>

#include <iostream>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include "cublas_v2.h"
#include "cutlass/gemm/device/gemm.h"

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
#define CUBLAS_CHECK(e)                                                \
  {                                                                    \
    cublasStatus_t err = (e);                                          \
    if (err != CUBLAS_STATUS_SUCCESS) {                                \
      fprintf(stderr, "CUBLAS error: %s, line %d, %s: %i\n", __FILE__, \
              __LINE__, #e, err);                                      \
      exit(EXIT_FAILURE);                                              \
    }                                                                  \
  }
#define CUTLASS_CHECK(e)                                                \
  {                                                                     \
    cutlass::Status err = (e);                                          \
    if (err != cutlass::Status::kSuccess) {                             \
      fprintf(stderr, "CUTLASS error: %s, line %d, %s: %i\n", __FILE__, \
              __LINE__, #e, (int)err);                                  \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
  }

namespace {
namespace detail {
struct pair_hash {
  std::size_t operator()(const std::pair<int, int> &p) const {
    return p.first * 10000 + p.second;
  }
};
}  // namespace detail

// this version is using cuda graphs and caches the graphs
template <typename Gemm, typename Real>
void run_group_graph(Gemm &&gemm, int m, int *n, int *k, Real alpha,
                     const Real *A, int lda, int *offsetsA, const Real *B,
                     int ldb, int *offsetsB, Real beta, Real *C, int ldc,
                     int *offsetsC, int batchCount, int blas_id = -1) {
  // we store at most one graph per "m" (# fields) and "blas id"
  static std::unordered_map<std::pair<int, int>, cudaGraphExec_t,
                            detail::pair_hash>
      graphCache;

  // we also store A, B, and C and recreate the graph if they change
  static std::unordered_map<
      std::pair<int, int>, std::tuple<Real const *, Real const *, Real const *>,
      detail::pair_hash>
      ptrCache;

  auto key = std::make_pair(m, blas_id);

  auto ptrs = ptrCache.find(key);
  if (ptrs != ptrCache.end() &&
      (std::get<0>(ptrs->second) != A || std::get<1>(ptrs->second) != B ||
       std::get<2>(ptrs->second) != C)) {
    // the plan is cached, but the pointers are not correct. we remove and
    // delete the graph, but we keep the cublas handles, if this happens more
    // often, we should cache this...
    std::cout << "WARNING: POINTER CHANGE --> THIS MIGHT BE SLOW" << std::endl;
    CUDA_CHECK(cudaGraphExecDestroy(graphCache[key]));
    graphCache.erase(key);
    ptrCache.erase(key);
  }

  auto graph = graphCache.find(key);
  if (graph == graphCache.end()) {
    // this graph does not exist yet
    cudaStream_t stream;
    CUDA_CHECK(cudaStreamCreate(&stream));

    cudaGraph_t new_graph;
    cudaGraphCreate(&new_graph, 0);
    for (int i = 0; i < batchCount; ++i) {
      if (m == 0 || n[i] == 0 || k[i] == 0) continue;

      CUDA_CHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal));
      gemm(stream, m, n[i], k[i], alpha, A + offsetsA[i], lda, B + offsetsB[i],
           ldb, beta, C + offsetsC[i], ldc);
      cudaGraph_t my_graph;
      CUDA_CHECK(cudaStreamEndCapture(stream, &my_graph));
      cudaGraphNode_t my_node;
      CUDA_CHECK(cudaGraphAddChildGraphNode(&my_node, new_graph, nullptr, 0,
                                            my_graph));
    }
    cudaGraphExec_t instance;
    CUDA_CHECK(cudaGraphInstantiate(&instance, new_graph, NULL, NULL, 0));
    CUDA_CHECK(cudaStreamDestroy(stream));
    CUDA_CHECK(cudaGraphDestroy(new_graph));

    graphCache.insert({key, instance});
    ptrCache.insert({key, std::make_tuple(A, B, C)});
  }

  CUDA_CHECK(cudaGraphLaunch(graphCache.at(key), 0));
}

// stupid simple gemm calls
template <typename Gemm, typename Real>
void run_group(Gemm &&gemm, int m, int *n, int *k, Real alpha, const Real *A,
               int lda, int *offsetsA, const Real *B, int ldb, int *offsetsB,
               Real beta, Real *C, int ldc, int *offsetsC, int batchCount,
               int = -1) {
  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0) continue;
    gemm(0, m, n[i], k[i], alpha, A + offsetsA[i], lda, B + offsetsB[i], ldb,
         beta, C + offsetsC[i], ldc);
  }
}

template <typename CutlassGemm>
CutlassGemm &get_cutlass_handle() {
  static auto handle = std::make_unique<CutlassGemm>();
  return *handle;
}

namespace detail {
template <cublasOperation_t TransA, cublasOperation_t TransB>
class cutlass_sgemm_grouped {
#if 0
  // we will enable this later (this ifdefs did not work, so I am going to enable this properly ltaer)
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
                         cutlass::layout::RowMajor>,  //
      float,
      std::conditional_t<TransB == CUBLAS_OP_N, cutlass::layout::ColumnMajor,
                         cutlass::layout::RowMajor>,  //
      float, cutlass::layout::ColumnMajor,            //
      float,                                          //
      OperatorClass, cutlass::arch::Sm80,             //
      ThreadblockShape, WarpShape, InstructionShape,  //
      cutlass::epilogue::thread::LinearCombination<   //
          float,                                      //
          128 / cutlass::sizeof_bits<float>::value,
          float,                                                     //
          float                                                      //
          >,                                                         //
      cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,  //
      3,                                                             //
      AlignmentA,                                                    //
      AlignmentB,                                                    //
      true,                                                          //
      MyOp                                                           //
      >;
  static constexpr int sz_align = 8;
#else
  // this was verified using Volta and uses FP32
  static constexpr int AlignmentA = 1;
  static constexpr int AlignmentB = 1;
  using ThreadblockShape = cutlass::gemm::GemmShape<128, 128, 8>;
  using WarpShape = cutlass::gemm::GemmShape<32, 32, 8>;
  using InstructionShape = cutlass::gemm::GemmShape<1, 1, 1>;
  using OperatorClass = cutlass::arch::OpClassSimt;
  using MyOp = cutlass::arch::OpMultiplyAdd;

  using Gemm = cutlass::gemm::device::Gemm<
      float,  //
      std::conditional_t<TransA == CUBLAS_OP_N, cutlass::layout::ColumnMajor,
                         cutlass::layout::RowMajor>,  //
      float,                                          //
      std::conditional_t<TransB == CUBLAS_OP_N, cutlass::layout::ColumnMajor,
                         cutlass::layout::RowMajor>,                 //
      float, cutlass::layout::ColumnMajor,                           //
      float,                                                         //
      OperatorClass, cutlass::arch::Sm70,                            //
      ThreadblockShape, WarpShape, InstructionShape,                 //
      cutlass::epilogue::thread::LinearCombination<                  //
          float,                                                     //
          1,                                                         //
          float,                                                     //
          float                                                      //
          >,                                                         //
      cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,  //
      2,                                                             //
      AlignmentA,                                                    //
      AlignmentB,                                                    //
      true,                                                          //
      MyOp                                                           //
      >;
  static constexpr int sz_align = 1;
#endif

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

}  // namespace detail
template <cublasOperation_t TransA, cublasOperation_t TransB>
void cutlass_sgemm_wrapper_grouped_op(int blas_id, int m, int *n, int *k,
                                      float alpha, const float *A, int lda,
                                      int *offsetsA, const float *B, int ldb,
                                      int *offsetsB, float beta, float *C,
                                      int ldc, int *offsetsC, int batchCount) {
  using namespace detail;
  run_group_graph(cutlass_sgemm_grouped<TransA, TransB>(), m, n, k, alpha, A,
                  lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
                  batchCount, blas_id);
}
void cutlass_sgemm_wrapper_grouped(int blas_id, cublasOperation_t transa,
                                   cublasOperation_t transb, int m, int *n,
                                   int *k, float alpha, const float *A, int lda,
                                   int *offsetsA, const float *B, int ldb,
                                   int *offsetsB, float beta, float *C, int ldc,
                                   int *offsetsC, int batchCount) {
  if (transa == CUBLAS_OP_N && transb == CUBLAS_OP_N)
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_N, CUBLAS_OP_N>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount);
  else if (transa == CUBLAS_OP_N && transb == CUBLAS_OP_T)
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_N, CUBLAS_OP_T>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount);
  else if (transa == CUBLAS_OP_T && transb == CUBLAS_OP_N)
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_T, CUBLAS_OP_N>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount);
  else if (transa == CUBLAS_OP_T && transb == CUBLAS_OP_T)
    cutlass_sgemm_wrapper_grouped_op<CUBLAS_OP_T, CUBLAS_OP_T>(
        blas_id, m, n, k, alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C,
        ldc, offsetsC, batchCount);
  else
    assert(false);
}

namespace detail {
cublasHandle_t get_cublas_handle() {
  static cublasHandle_t handle;
  if (!handle) CUBLAS_CHECK(cublasCreate(&handle));
  return handle;
}
template <typename Real>
struct cublas_gemm_grouped {
 public:
  cublas_gemm_grouped(cublasOperation_t transa, cublasOperation_t transb)
      : transa_(transa), transb_(transb) {
    // we need to get the cublas handle here, otherwise this could be created
    // during graph capturing
    get_cublas_handle();
  };
  void operator()(cudaStream_t stream, int m, int n, int k, Real alpha,
                  const Real *A, int lda, const Real *B, int ldb, Real beta,
                  Real *C, int ldc) const {
    cublasHandle_t handle = get_cublas_handle();
    CUBLAS_CHECK(cublasSetStream(handle, stream));

    if constexpr (std::is_same<Real, float>::value)
      CUBLAS_CHECK(cublasSgemm(handle, transa_, transb_, m, n, k, &alpha, A,
                               lda, B, ldb, &beta, C, ldc));
    if constexpr (std::is_same<Real, double>::value)
      CUBLAS_CHECK(cublasDgemm(handle, transa_, transb_, m, n, k, &alpha, A,
                               lda, B, ldb, &beta, C, ldc));
  }

 private:
  cublasOperation_t transa_, transb_;
};
}  // namespace detail
void cublas_sgemm_wrapper_grouped(int blas_id, cublasOperation_t transa,
                                  cublasOperation_t transb, int m, int *n,
                                  int *k, float alpha, const float *A, int lda,
                                  int *offsetsA, const float *B, int ldb,
                                  int *offsetsB, float beta, float *C, int ldc,
                                  int *offsetsC, int batchCount) {
  using namespace detail;
  run_group_graph(cublas_gemm_grouped<float>(transa, transb), m, n, k, alpha, A,
                  lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
                  batchCount, blas_id);
}
void cublas_dgemm_wrapper_grouped(int blas_id, cublasOperation_t transa,
                                  cublasOperation_t transb, int m, int *n,
                                  int *k, double alpha, const double *A,
                                  int lda, int *offsetsA, const double *B,
                                  int ldb, int *offsetsB, double beta,
                                  double *C, int ldc, int *offsetsC,
                                  int batchCount) {
  using namespace detail;
  run_group_graph(cublas_gemm_grouped<double>(transa, transb), m, n, k, alpha,
                  A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
                  batchCount, blas_id);
}
}  // namespace

extern "C" {
void cublas_dgemm_wrapper(cublasOperation_t transa, cublasOperation_t transb,
                          int m, int n, int k, double alpha, const double *A,
                          int lda, int tda, const double *B, int ldb, int tdb,
                          double beta, double *C, int ldc, int tdc,
                          int batchCount) {
  static cublasHandle_t handle = nullptr;
  if (!handle) CUBLAS_CHECK(cublasCreate(&handle));

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
  if (!handle) CUBLAS_CHECK(cublasCreate(&handle));

  CUBLAS_CHECK(cublasSgemmStridedBatched(handle, transa, transb, m, n, k,
                                         &alpha, A, lda, tda, B, ldb, tdb,
                                         &beta, C, ldc, tdc, batchCount));
}

void blas_sgemm_wrapper_grouped(int blas_id, cublasOperation_t transa,
                                cublasOperation_t transb, int m, int *n, int *k,
                                float alpha, const float *A, int lda,
                                int *offsetsA, const float *B, int ldb,
                                int *offsetsB, float beta, float *C, int ldc,
                                int *offsetsC, int batchCount) {
  if (use_cutlass)
    cutlass_sgemm_wrapper_grouped(blas_id, transa, transb, m, n, k, alpha, A,
                                  lda, offsetsA, B, ldb, offsetsB, beta, C, ldc,
                                  offsetsC, batchCount);
  else
    cublas_sgemm_wrapper_grouped(blas_id, transa, transb, m, n, k, alpha, A, lda,
                                 offsetsA, B, ldb, offsetsB, beta, C, ldc,
                                 offsetsC, batchCount);
}
void blas_dgemm_wrapper_grouped(int blas_id, cublasOperation_t transa,
                                cublasOperation_t transb, int m, int *n, int *k,
                                double alpha, const double *A, int lda,
                                int *offsetsA, const double *B, int ldb,
                                int *offsetsB, double beta, double *C, int ldc,
                                int *offsetsC, int batchCount) {
  cublas_dgemm_wrapper_grouped(blas_id, transa, transb, m, n, k, alpha, A, lda, offsetsA,
                               B, ldb, offsetsB, beta, C, ldc, offsetsC,
                               batchCount);
}
}
