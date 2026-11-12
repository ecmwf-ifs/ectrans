// (C) Copyright 2026- ECMWF.
// Licensed under the Apache License, Version 2.0.
#pragma once

// Include hicblas.h first. Scheduling is shared by dense GEMMs and butterfly phases.
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <type_traits>
#include <vector>

namespace gemm_execution {

inline int setting(const char *name, int fallback) {
  const char *value = std::getenv(name);
  if (!value) return fallback;
  char *end;
  errno = 0;
  long n = std::strtol(value, &end, 10);
  if (errno || end == value || *end || n < 0 || n > INT_MAX) {
    std::fprintf(stderr, "%s must be a nonnegative integer\n", name);
    std::abort();
  }
  return int(n);
}

struct Options {
  bool graphs;
  int branches; // 0: independent work; positive: bounded stream chains.
};

inline Options options(bool butterfly, int default_branches) {
#ifdef USE_GRAPHS_GEMM
  constexpr bool default_graphs = true;
#else
  constexpr bool default_graphs = false;
#endif
  int graphs = butterfly
      ? setting("ECTRANS_BUTTERFLY_GRAPHS", default_graphs)
      : setting("ECTRANS_DENSE_GRAPHS", default_graphs);
  if (graphs > 1) {
    std::fprintf(stderr, "ECTRANS_%s_GRAPHS must be 0 or 1\n", butterfly ? "BUTTERFLY" : "DENSE");
    std::abort();
  }
  return {graphs != 0, butterfly
      ? setting("ECTRANS_BUTTERFLY_NBRANCHES", default_branches)
      : setting("ECTRANS_DENSE_NBRANCHES", default_branches)};
}

// One scheduler per workload/precision. Resources persist like the graph caches;
// reserve/prepare must run before capture, never while constructing graph nodes.
class Scheduler {
  std::vector<hipStream_t> streams_;
  std::vector<hipEvent_t> joins_;
  hipEvent_t fork_;
  hipStream_t capture_ = nullptr;
  bool initialized_ = false;

public:
  const Options config;
  explicit Scheduler(Options config) : config(config) {}
  Scheduler(const Scheduler &) = delete;
  Scheduler &operator=(const Scheduler &) = delete;

  int width(int tasks) const { return config.branches ? std::min(tasks, config.branches) : tasks; }
  hipStream_t capture_stream() const { return capture_; }

  void reserve(int tasks) {
    if (!initialized_) {
      HIC_CHECK(hipStreamCreateWithFlags(&capture_, hipStreamNonBlocking));
      HIC_CHECK(hipEventCreateWithFlags(&fork_, hipEventDisableTiming));
      initialized_ = true;
    }
    while (int(streams_.size()) < width(tasks)) {
      hipStream_t stream;
      hipEvent_t event;
      HIC_CHECK(hipStreamCreateWithFlags(&stream, hipStreamNonBlocking));
      HIC_CHECK(hipEventCreateWithFlags(&event, hipEventDisableTiming));
      streams_.push_back(stream);
      joins_.push_back(event);
    }
  }

  template <typename Prepare> void prepare(int tasks, hipStream_t caller, Prepare prepare) {
    reserve(tasks);
    prepare(caller);
    prepare(capture_);
    for (int i = 0; i < width(tasks); ++i) prepare(streams_[i]);
  }

  // The same fork/join works for direct execution and stream capture.
  template <typename Launch> void run(int tasks, hipStream_t stream, Launch launch) const {
    const int count = width(tasks);
    if (count <= 1) {
      for (int i = 0; i < tasks; ++i) launch(i, stream);
      return;
    }
    HIC_CHECK(hipEventRecord(fork_, stream));
    for (int i = 0; i < count; ++i) HIC_CHECK(hipStreamWaitEvent(streams_[i], fork_, 0));
    for (int i = 0; i < tasks; ++i) launch(i, streams_[i % count]);
    for (int i = 0; i < count; ++i) {
      HIC_CHECK(hipEventRecord(joins_[i], streams_[i]));
      HIC_CHECK(hipStreamWaitEvent(stream, joins_[i], 0));
    }
  }

  template <typename Launch>
  hipGraphNode_t append_serial(hipGraph_t graph, const std::vector<hipGraphNode_t> &parents,
                              hipStream_t stream, Launch launch) const {
    HIC_CHECK(hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal));
    launch(stream);
    hipGraph_t child;
    HIC_CHECK(hipStreamEndCapture(stream, &child));
    hipGraphNode_t node;
    HIC_CHECK(hipGraphAddChildGraphNode(&node, graph, parents.data(), parents.size(), child));
    HIC_CHECK(hipGraphDestroy(child));
    return node;
  }

  // A phase joins all children before its successor, without ordering siblings.
  template <typename Launch>
  void append(hipGraph_t graph, std::vector<hipGraphNode_t> &parents,
              int tasks, Launch launch) const {
    if (!tasks) return;
    std::vector<hipGraphNode_t> children;
    for (int i = 0; i < tasks; ++i)
      children.push_back(append_serial(graph, parents, streams_[i],
                                      [&](hipStream_t stream) { launch(i, stream); }));
    parents = std::move(children);
  }

  template <typename Launch> hipGraph_t graph(int tasks, Launch launch) const {
    hipGraph_t graph;
    if (config.branches == 0) {
      HIC_CHECK(hipGraphCreate(&graph, 0));
      std::vector<hipGraphNode_t> parents;
      append(graph, parents, tasks, launch);
    } else {
      HIC_CHECK(hipStreamBeginCapture(capture_, hipStreamCaptureModeGlobal));
      run(tasks, capture_, launch);
      HIC_CHECK(hipStreamEndCapture(capture_, &graph));
    }
    return graph;
  }
};

// Host enqueueing is serial. CUDA shares a handle; HIP needs one per stream.
// https://rocm.docs.amd.com/projects/rocBLAS/en/docs-7.2.0/how-to/Programmers_Guide.html
class BlasContext {
#ifdef CUDAGPU
  hipblasHandle_t shared_ = nullptr;
  std::map<hipStream_t, void *> workspaces_;
#else
  std::map<hipStream_t, hipblasHandle_t> handles_;
#endif
  std::map<hipStream_t, unsigned> warmed_;
public:
  hipblasHandle_t handle(hipStream_t stream) {
#ifdef CUDAGPU
    if (!shared_) HICBLAS_CHECK(hipblasCreate(&shared_));
    hipblasHandle_t handle = shared_;
#else
    auto &handle = handles_[stream];
    if (!handle) HICBLAS_CHECK(hipblasCreate(&handle));
#endif
    HICBLAS_CHECK(hipblasSetStream(handle, stream));
#ifdef CUDAGPU
    auto it = workspaces_.find(stream);
    if (it != workspaces_.end()) HICBLAS_CHECK(cublasSetWorkspace(handle, it->second, 32 * 1024 * 1024));
#endif
    return handle;
  }

  template <typename Real> void prepare(hipStream_t stream, bool children) {
#ifdef CUDAGPU
    // Child graphs cannot contain cuBLAS allocation nodes. Give each concurrent
    // stream its own workspace; 32 MiB is NVIDIA's Hopper recommendation.
    // https://docs.nvidia.com/cuda/archive/13.0.0/cublas/index.html#cuda-graphs-support
    if (children && !workspaces_.count(stream)) {
      void *workspace;
      HIC_CHECK(hipMalloc(&workspace, 32 * 1024 * 1024));
      workspaces_[stream] = workspace;
    }
#endif
    const unsigned bit = std::is_same<Real, float>::value ? 1 : 2;
    if (warmed_[stream] & bit) return;
    Real *scratch;
    HIC_CHECK(hipMalloc(&scratch, 3 * sizeof(Real)));
    HIC_CHECK(hipMemset(scratch, 0, 3 * sizeof(Real)));
    const Real one = 1, zero = 0;
    auto h = handle(stream);
    if constexpr (std::is_same<Real, float>::value) {
      HICBLAS_CHECK(hipblasSgemm(h, HIPBLAS_OP_N, HIPBLAS_OP_N, 1, 1, 1,
                               &one, scratch, 1, scratch + 1, 1, &zero, scratch + 2, 1));
    } else {
      HICBLAS_CHECK(hipblasDgemm(h, HIPBLAS_OP_N, HIPBLAS_OP_N, 1, 1, 1,
                               &one, scratch, 1, scratch + 1, 1, &zero, scratch + 2, 1));
    }
    HIC_CHECK(hipStreamSynchronize(stream));
    HIC_CHECK(hipFree(scratch));
    warmed_[stream] |= bit;
  }
};
inline BlasContext &blas() { static BlasContext context; return context; }

} // namespace gemm_execution
