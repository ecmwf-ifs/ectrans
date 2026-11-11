// (C) Copyright 2026- NVIDIA.
// Licensed under the Apache License, Version 2.0.
#pragma once

#include <cublasLt.h>
#include <cuda_runtime.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <tuple>

namespace legendre_cublaslt {

inline void check(cublasStatus_t status) {
  if (status != CUBLAS_STATUS_SUCCESS) {
    std::fprintf(stderr, "Legendre cuBLASLt: status %d\n", int(status));
    std::abort();
  }
}

inline void check(cudaError_t status) {
  if (status != cudaSuccess) {
    std::fprintf(stderr, "Legendre cuBLASLt: %s\n", cudaGetErrorString(status));
    std::abort();
  }
}

// Read once: changing the backend after graph capture is unsupported.
inline bool enabled() {
  static const bool value = [] {
    const char *backend = std::getenv("ECTRANS_BLAS_IMPLEMENTATION");
    return backend && std::strcmp(backend, "cublaslt_fp32") == 0;
  }();
  return value;
}

struct Handle {
  cublasLtHandle_t value = nullptr;
  Handle() { check(cublasLtCreate(&value)); }
  ~Handle() { cublasLtDestroy(value); }
  Handle(const Handle &) = delete;
  Handle &operator=(const Handle &) = delete;
};

inline cublasLtHandle_t handle() {
  static Handle instance;
  return instance.value;
}

inline uint32_t alignment(const void *pointer) {
  uint32_t bytes = 256;
  while (reinterpret_cast<uintptr_t>(pointer) % bytes != 0)
    bytes /= 2;
  return bytes;
}

// Immutable layout/algorithm identity. Each MPI process owns one GPU.
struct Shape {
  bool transpose_a, transpose_b;
  int m, n, k, lda, ldb, ldc, batch_count;
  bool pointer_array;
  uint32_t alignment_a, alignment_b, alignment_c;

  auto key() const {
    return std::tie(transpose_a, transpose_b, m, n, k, lda, ldb, ldc,
                    batch_count, pointer_array, alignment_a, alignment_b,
                    alignment_c);
  }
  bool operator<(const Shape &other) const { return key() < other.key(); }
};

inline Shape shape(bool transpose_a, bool transpose_b, int m, int n, int k,
                   int lda, int ldb, int ldc, int batch_count, bool pointer_array,
                   const void *a, const void *b, const void *c) {
  // A device pointer array's own alignment says nothing about its entries.
  // FP32 entries are guaranteed only four-byte alignment.
  return {transpose_a, transpose_b, m, n, k, lda, ldb, ldc, batch_count,
          pointer_array, pointer_array ? 4 : alignment(a),
          pointer_array ? 4 : alignment(b), pointer_array ? 4 : alignment(c)};
}

struct Plan {
  cublasLtMatmulDesc_t operation = nullptr;
  cublasLtMatrixLayout_t a = nullptr, b = nullptr, c = nullptr;
  cublasLtMatmulAlgo_t algorithm;

  explicit Plan(const Shape &s) {
    check(cublasLtMatmulDescCreate(&operation, CUBLAS_COMPUTE_32F_PEDANTIC,
                                  CUDA_R_32F));
    cublasOperation_t ta = s.transpose_a ? CUBLAS_OP_T : CUBLAS_OP_N;
    cublasOperation_t tb = s.transpose_b ? CUBLAS_OP_T : CUBLAS_OP_N;
    check(cublasLtMatmulDescSetAttribute(operation, CUBLASLT_MATMUL_DESC_TRANSA,
                                        &ta, sizeof(ta)));
    check(cublasLtMatmulDescSetAttribute(operation, CUBLASLT_MATMUL_DESC_TRANSB,
                                        &tb, sizeof(tb)));
    check(cublasLtMatrixLayoutCreate(&a, CUDA_R_32F, s.transpose_a ? s.k : s.m,
                                    s.transpose_a ? s.m : s.k, s.lda));
    check(cublasLtMatrixLayoutCreate(&b, CUDA_R_32F, s.transpose_b ? s.n : s.k,
                                    s.transpose_b ? s.k : s.n, s.ldb));
    check(cublasLtMatrixLayoutCreate(&c, CUDA_R_32F, s.m, s.n, s.ldc));
    if (s.pointer_array) {
      int mode = CUBLASLT_BATCH_MODE_POINTER_ARRAY;
      for (auto layout : {a, b, c}) {
        check(cublasLtMatrixLayoutSetAttribute(
            layout, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &s.batch_count,
            sizeof(s.batch_count)));
        check(cublasLtMatrixLayoutSetAttribute(
            layout, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode)));
      }
    }

    // Keep the default heuristic choice. No workspace is shared by concurrent
    // child graphs or FLT branch streams; all GEMMs have zero workspace.
    cublasLtMatmulPreference_t preference;
    check(cublasLtMatmulPreferenceCreate(&preference));
    size_t workspace_bytes = 0;
    check(cublasLtMatmulPreferenceSetAttribute(
        preference, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &workspace_bytes,
        sizeof(workspace_bytes)));
    const uint32_t alignments[] = {s.alignment_a, s.alignment_b,
                                  s.alignment_c, s.alignment_c};
    int i = 0;
    for (auto attribute : {CUBLASLT_MATMUL_PREF_MIN_ALIGNMENT_A_BYTES,
                           CUBLASLT_MATMUL_PREF_MIN_ALIGNMENT_B_BYTES,
                           CUBLASLT_MATMUL_PREF_MIN_ALIGNMENT_C_BYTES,
                           CUBLASLT_MATMUL_PREF_MIN_ALIGNMENT_D_BYTES}) {
      check(cublasLtMatmulPreferenceSetAttribute(
          preference, attribute, &alignments[i++], sizeof(uint32_t)));
    }
    cublasLtMatmulHeuristicResult_t candidates[8];
    int found = 0;
    check(cublasLtMatmulAlgoGetHeuristic(handle(), operation, a, b, c, c,
                                        preference, 8, candidates, &found));
    check(cublasLtMatmulPreferenceDestroy(preference));
    for (int j = 0; j < found; ++j) {
      if (candidates[j].state == CUBLAS_STATUS_SUCCESS &&
          candidates[j].workspaceSize == 0) {
        algorithm = candidates[j].algo;
        return;
      }
    }
    std::fprintf(stderr, "Legendre cuBLASLt: no FP32 zero-workspace algorithm "
                         "(m=%d n=%d k=%d lda=%d ldb=%d ldc=%d batch=%d)\n",
                 s.m, s.n, s.k, s.lda, s.ldb, s.ldc, s.batch_count);
    std::abort();
  }

  ~Plan() {
    if (operation) cublasLtMatmulDescDestroy(operation);
    if (a) cublasLtMatrixLayoutDestroy(a);
    if (b) cublasLtMatrixLayoutDestroy(b);
    if (c) cublasLtMatrixLayoutDestroy(c);
  }
  Plan(const Plan &) = delete;
  Plan &operator=(const Plan &) = delete;
};

inline Plan &plan(const Shape &s, cudaStream_t stream = nullptr) {
  static std::map<Shape, std::unique_ptr<Plan>> cache;
  auto found = cache.find(s);
  if (found == cache.end()) {
    cudaStreamCaptureStatus status;
    check(cudaStreamIsCapturing(stream, &status));
    if (status != cudaStreamCaptureStatusNone) {
      std::fprintf(stderr, "Legendre cuBLASLt: prepare plans before capture\n");
      std::abort();
    }
    found = cache.emplace(s, std::make_unique<Plan>(s)).first;
  }
  return *found->second;
}

inline void prepare(bool transpose_a, bool transpose_b, int m, int n, int k,
                    int lda, int ldb, int ldc, int batch_count, bool pointer_array,
                    const void *a, const void *b, const void *c) {
  if (m <= 0 || n <= 0 || k <= 0 || batch_count <= 0) return;
  (void)plan(shape(transpose_a, transpose_b, m, n, k, lda, ldb, ldc, batch_count,
                   pointer_array, a, b, c));
}

// Single and pointer-array GEMMs share the same descriptor cache and launch.
inline void gemm(cudaStream_t stream, bool transpose_a, bool transpose_b,
                 int m, int n, int k, float alpha, const void *a, int lda,
                 const void *b, int ldb, float beta, void *c, int ldc,
                 int batch_count = 1, bool pointer_array = false) {
  if (m <= 0 || n <= 0 || k <= 0 || batch_count <= 0) return;
  auto &p = plan(shape(transpose_a, transpose_b, m, n, k, lda, ldb, ldc,
                       batch_count, pointer_array, a, b, c), stream);
  check(cublasLtMatmul(handle(), p.operation, &alpha, a, p.a, b, p.b, &beta,
                       c, p.c, c, p.c, &p.algorithm, nullptr, 0, stream));
}

} // namespace legendre_cublaslt
