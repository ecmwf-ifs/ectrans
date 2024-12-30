// (C) Copyright 2000- ECMWF.
// (C) Copyright 2024- NVIDIA.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include "hicblas.h"
#ifdef USE_CUTLASS
#include "cutlass/gemm/device/gemm.h"
#endif

#include "growing_allocator.h"

bool hip_alreadyAllocated_sgemm = false;
bool hip_alreadyAllocated_sgemm_handle = false;

bool hip_alreadyAllocated_dsgemm = false;
bool hip_alreadyAllocated_dgemm_handle = false;

hipblasHandle_t handle_hip_sgemm;
hipblasHandle_t handle_hip_dgemm;

namespace {
struct cache_key {
  int resol_id;
  int m;
  int blas_id;

  bool operator==(const cache_key &other) const {
    return resol_id == other.resol_id && m == other.m &&
           blas_id == other.blas_id;
  }
  cache_key(int resol_id_, int m_, int blas_id_)
      : resol_id(resol_id_), m(m_), blas_id(blas_id_) {}
};
} // namespace
template <> struct std::hash<cache_key> {
  std::size_t operator()(const cache_key &k) const {
    return k.blas_id * 1000000 + k.resol_id * 10000 + k.m;
  }
};

namespace {
template <typename Gemm> auto &get_graph_cache() {
  // we store at most one graph per "m" (# fields) and "blas id" and resolution
  static std::unordered_map<cache_key, std::shared_ptr<hipGraphExec_t>>
      graphCache;
  return graphCache;
}
template <typename Gemm> auto &get_ptr_cache() {
  using real_t = typename Gemm::real_type;
  static std::unordered_map<
      cache_key, std::tuple<real_t const *, real_t const *, real_t const *>>
      ptrCache;
  return ptrCache;
}

template <typename Gemm> void free_gemm_graph_cache(float *, size_t) {
  get_graph_cache<Gemm>().clear();
  get_ptr_cache<Gemm>().clear();
}
template <typename Cache>
void erase_resol_from_cache(Cache &cache, int resol_id) {
  // Note that in C++20 this could also be std::erase_if
  int erased = 0;
  for (auto it = cache.begin(); it != cache.end();) {
    if (it->first.resol_id == resol_id) {
      it = cache.erase(it);
      ++erased;
    } else
      ++it;
  }
}
template <typename Gemm> void erase_from_caches(int resol_id) {
  erase_resol_from_cache(get_graph_cache<Gemm>(), resol_id);
  erase_resol_from_cache(get_ptr_cache<Gemm>(), resol_id);
}

// this version is using graphs and caches the graphs
template <typename Gemm, typename Real>
void run_group_graph(Gemm &&gemm, int resol_id, int m, const int *n,
                     const int *k, Real alpha, const Real *A, int lda,
                     const int64_t *offsetsA, const Real *B, const int *ldb,
                     const int64_t *offsetsB, Real beta, Real *C, int ldc,
                     const int64_t *offsetsC, int batchCount,
                     hipStream_t stream, int blas_id, void *growing_allocator) {
  growing_allocator_register_free_c(growing_allocator,
                                    free_gemm_graph_cache<Gemm>);

  // we store at most one graph per "m" (# fields) and "blas id"
  auto &graphCache = get_graph_cache<Gemm>();

  // we also store A, B, and C and recreate the graph if they change
  auto &ptrCache = get_ptr_cache<Gemm>();

  auto key = cache_key{resol_id, m, blas_id};

  auto ptrs = ptrCache.find(key);
  if (ptrs != ptrCache.end() &&
      (std::get<0>(ptrs->second) != A || std::get<1>(ptrs->second) != B ||
       std::get<2>(ptrs->second) != C)) {
    // the plan is cached, but the pointers are not correct. we remove and
    // delete the graph, but we keep the hipblas handles, if this happens more
    // often, we should cache this...
    std::cout
        << "WARNING GEMM: POINTER CHANGE - Graph recreation might be slow."
        << std::endl;
    std::cout << "We have an entry with key {m=" << m << ", blas_id=" << blas_id
              << ", resol_id=" << resol_id << "}\n";
    std::cout << "Pointers: " << std::get<0>(ptrs->second) << ", "
              << std::get<1>(ptrs->second) << ", " << std::get<2>(ptrs->second)
              << " vs. " << A << ", " << B << ", " << C << std::endl;
    graphCache.erase(key);
    ptrCache.erase(key);
  }

  auto graph = graphCache.find(key);
  if (graph == graphCache.end()) {
    // this graph does not exist yet
    hipStream_t stream;
    HIC_CHECK(hipStreamCreate(&stream));

    hipGraph_t new_graph;
    hipGraphCreate(&new_graph, 0);
    for (int i = 0; i < batchCount; ++i) {
      if (m == 0 || n[i] == 0 || k[i] == 0)
        continue;

      HIC_CHECK(hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal));
      gemm(stream, m, n[i], k[i], alpha, A + offsetsA[i], lda, B + offsetsB[i],
           ldb[i], beta, C + offsetsC[i], ldc);
      hipGraph_t my_graph;
      HIC_CHECK(hipStreamEndCapture(stream, &my_graph));
      hipGraphNode_t my_node;
      HIC_CHECK(
          hipGraphAddChildGraphNode(&my_node, new_graph, nullptr, 0, my_graph));
    }
    hipGraphExec_t instance;
    HIC_CHECK(hipGraphInstantiate(&instance, new_graph, NULL, NULL, 0));
    HIC_CHECK(hipStreamDestroy(stream));
    HIC_CHECK(hipGraphDestroy(new_graph));

    graphCache.insert({key, std::shared_ptr<hipGraphExec_t>(
                                new hipGraphExec_t{instance}, [](auto ptr) {
                                  HIC_CHECK(hipGraphExecDestroy(*ptr));
                                  delete ptr;
                                })});
    ptrCache.insert({key, std::make_tuple(A, B, C)});
  }

  HIC_CHECK(hipGraphLaunch(*graphCache.at(key), stream));
}

// stupid simple gemm calls
template <typename Gemm, typename Real>
void run_group(Gemm &&gemm, int resol_id, int m, const int *n, const int *k,
               Real alpha, const Real *A, int lda, const int64_t *offsetsA,
               const Real *B, const int *ldb, const int64_t *offsetsB,
               Real beta, Real *C, int ldc, const int64_t *offsetsC,
               int batchCount, hipStream_t stream, int = -1) {
  for (int i = 0; i < batchCount; ++i) {
    if (m == 0 || n[i] == 0 || k[i] == 0)
      continue;
    gemm(stream, m, n[i], k[i], alpha, A + offsetsA[i], lda, B + offsetsB[i],
         ldb[i], beta, C + offsetsC[i], ldc);
  }
}

#ifdef USE_CUTLASS
#include "hicblas_cutlass.cuda.h"
#endif

hipblasHandle_t get_hipblas_handle() {
  static hipblasHandle_t handle;
  if (!handle)
    HICBLAS_CHECK(hipblasCreate(&handle));
  return handle;
}
template <typename Real> struct hipblas_gemm_grouped {
public:
  using real_type = Real;
  hipblas_gemm_grouped(hipblasOperation_t transa, hipblasOperation_t transb)
      : transa_(transa), transb_(transb) {
    // we need to get the hipblas handle here, otherwise this could be created
    // during graph capturing
    get_hipblas_handle();
  };
  void operator()(hipStream_t stream, int m, int n, int k, Real alpha,
                  const Real *A, int lda, const Real *B, int ldb, Real beta,
                  Real *C, int ldc) const {
    hipblasHandle_t handle = get_hipblas_handle();
    HICBLAS_CHECK(hipblasSetStream(handle, stream));

    if constexpr (std::is_same<Real, float>::value)
      HICBLAS_CHECK(hipblasSgemm(handle, transa_, transb_, m, n, k, &alpha, A,
                                 lda, B, ldb, &beta, C, ldc));
    if constexpr (std::is_same<Real, double>::value)
      HICBLAS_CHECK(hipblasDgemm(handle, transa_, transb_, m, n, k, &alpha, A,
                                 lda, B, ldb, &beta, C, ldc));
  }

private:
  hipblasOperation_t transa_, transb_;
};

#ifndef USE_CUTLASS

void hipblas_sgemm_wrapper_grouped(
    int resol_id, int blas_id, char transa, char transb, int m, const int *n,
    const int *k, float alpha, const float *A, int lda, const int64_t *offsetsA,
    const float *B, const int *ldb, const int64_t *offsetsB, float beta,
    float *C, int ldc, const int64_t *offsetsC, int batchCount,
    hipStream_t stream, void *growing_allocator) {

  hipblasOperation_t op_t1 = HIPBLAS_OP_N, op_t2 = HIPBLAS_OP_N;
  if (transa == 'T' || transa == 't')
    op_t1 = HIPBLAS_OP_T;
  if (transb == 'T' || transb == 't')
    op_t2 = HIPBLAS_OP_T;

#ifdef USE_GRAPHS_GEMM
  run_group_graph(hipblas_gemm_grouped<float>(op_t1, op_t2), resol_id, m, n, k,
                  alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc,
                  offsetsC, batchCount, stream, blas_id, growing_allocator);
#else
  run_group(hipblas_gemm_grouped<float>(op_t1, op_t2), resol_id, m, n, k, alpha,
            A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
            batchCount, stream);
#endif
}

#endif

void hipblas_dgemm_wrapper_grouped(int resol_id, int blas_id, char transa,
                                   char transb, int m, const int *n,
                                   const int *k, double alpha, const double *A,
                                   int lda, const int64_t *offsetsA,
                                   double const *B, const int *ldb,
                                   const int64_t *offsetsB, double beta,
                                   double *C, int ldc, const int64_t *offsetsC,
                                   int batchCount, hipStream_t stream, void *) {

  hipblasOperation_t op_t1 = HIPBLAS_OP_N, op_t2 = HIPBLAS_OP_N;
  if (transa == 'T' || transa == 't')
    op_t1 = HIPBLAS_OP_T;
  if (transb == 'T' || transb == 't')
    op_t2 = HIPBLAS_OP_T;

  run_group(hipblas_gemm_grouped<double>(op_t1, op_t2), resol_id, m, n, k,
            alpha, A, lda, offsetsA, B, ldb, offsetsB, beta, C, ldc, offsetsC,
            batchCount, stream, blas_id);
}

} // namespace

extern "C" {
void hipblas_dgemm_wrapper(char transa, char transb, int m, int n, int k,
                           double alpha, const double *A, int lda, int tda,
                           const double *B, int ldb, int tdb, double beta,
                           double *C, int ldc, int tdc, int batchCount,
                           size_t stream, void *growing_allocator) {

  hipblasOperation_t op_t1 = HIPBLAS_OP_N, op_t2 = HIPBLAS_OP_N;

  if (transa == 'T' || transa == 't')
    op_t1 = HIPBLAS_OP_T;
  if (transb == 'T' || transb == 't')
    op_t2 = HIPBLAS_OP_T;

  if (!hip_alreadyAllocated_dgemm_handle) {
    HICBLAS_CHECK(hipblasCreate(&handle_hip_dgemm));
    hip_alreadyAllocated_dgemm_handle = true;
  }
  hipblasHandle_t handle = get_hipblas_handle();
  HICBLAS_CHECK(hipblasSetStream(handle, *(hipStream_t *)stream));

  HICBLAS_CHECK(hipblasDgemmStridedBatched(
      handle, op_t1, op_t2, m, n, k, &alpha, (const double *)A, lda, tda,
      (const double *)B, ldb, tdb, &beta, (double *)C, ldc, tdc, batchCount));
}

void hipblas_sgemm_wrapper(char transa, char transb, int m, int n, int k,
                           float alpha, const float *A, int lda, int tda,
                           const float *B, int ldb, int tdb, float beta,
                           float *C, int ldc, int tdc, int batchCount,
                           void *growing_allocator) {

  hipblasOperation_t op_t1 = HIPBLAS_OP_N, op_t2 = HIPBLAS_OP_N;

  if (transa == 'T' || transa == 't')
    op_t1 = HIPBLAS_OP_T;
  if (transb == 'T' || transb == 't')
    op_t2 = HIPBLAS_OP_T;

  if (!hip_alreadyAllocated_sgemm_handle) {
    HICBLAS_CHECK(hipblasCreate(&handle_hip_sgemm));
    hip_alreadyAllocated_sgemm_handle = true;
  }
  HICBLAS_CHECK(hipblasSgemmStridedBatched(
      handle_hip_sgemm, op_t1, op_t2, m, n, k, &alpha, (const float *)A, lda,
      tda, (const float *)B, ldb, tdb, &beta, (float *)C, ldc, tdc,
      batchCount));
}

void hipblas_sgemm_wrapper_grouped(
    int resol_id, int blas_id, char transa, char transb, int m, const int *n,
    const int *k, float alpha, const float *A, int lda, const int64_t *offsetsA,
    const float *B, const int *ldb, const int64_t *offsetsB, float beta,
    float *C, int ldc, const int64_t *offsetsC, int batchCount, size_t stream,
    void *growing_allocator) {
#ifdef USE_CUTLASS
  cutlass_sgemm_wrapper_grouped(resol_id, blas_id, transa, transb, m, n, k,
                                alpha, A, lda, offsetsA, B, ldb, offsetsB, beta,
                                C, ldc, offsetsC, batchCount,
                                *(hipStream_t *)stream, growing_allocator);
#else
  hipblas_sgemm_wrapper_grouped(resol_id, blas_id, transa, transb, m, n, k,
                                alpha, A, lda, offsetsA, B, ldb, offsetsB, beta,
                                C, ldc, offsetsC, batchCount,
                                *(hipStream_t *)stream, growing_allocator);
#endif
}

void hipblas_dgemm_wrapper_grouped(int resol_id, int blas_id, char transa,
                                   char transb, int m, const int *n,
                                   const int *k, double alpha, double const *A,
                                   int lda, const int64_t *offsetsA,
                                   double const *B, const int *ldb,
                                   const int64_t *offsetsB, double beta,
                                   double *C, int ldc, const int64_t *offsetsC,
                                   int batchCount, size_t stream,
                                   void *growing_allocator) {
  hipblas_dgemm_wrapper_grouped(resol_id, blas_id, transa, transb, m, n, k,
                                alpha, A, lda, offsetsA, B, ldb, offsetsB, beta,
                                C, ldc, offsetsC, batchCount,
                                *(hipStream_t *)stream, growing_allocator);
}

void clean_gemm(int resol_id) {
  erase_from_caches<hipblas_gemm_grouped<float>>(resol_id);
  erase_from_caches<hipblas_gemm_grouped<double>>(resol_id);
#ifdef USE_CUTLASS
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_3xtf32, CUBLAS_OP_T, CUBLAS_OP_T>>(resol_id);
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_3xtf32, CUBLAS_OP_N, CUBLAS_OP_T>>(resol_id);
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_3xtf32, CUBLAS_OP_T, CUBLAS_OP_N>>(resol_id);
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_3xtf32, CUBLAS_OP_N, CUBLAS_OP_N>>(resol_id);
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_fp32, CUBLAS_OP_T, CUBLAS_OP_T>>(resol_id);
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_fp32, CUBLAS_OP_N, CUBLAS_OP_T>>(resol_id);
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_fp32, CUBLAS_OP_T, CUBLAS_OP_N>>(resol_id);
  erase_from_caches<detail::cutlass_sgemm_grouped<
      detail::CutlassType::cutlass_fp32, CUBLAS_OP_N, CUBLAS_OP_N>>(resol_id);
#endif
}
}
