// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

// GPU version of MULT_BUTM from butterfly_alg_mod.F90, operating on the flattened butterfly
// structure (BUTTERFLY_STRUCT_FLAT) for all FLT zonal wavenumbers in one call.
//
// The input (A) and output (C) arrays are stored with the field index fastest, as for the other
// GPU Legendre GEMMs, i.e. element (row, field) of mode m is at A[a_offsets[m] + field + row * lda].
// The internal work arrays (beta and vec) are stored the same way, with leading dimension n_flds,
// i.e. they hold the transposes of the Fortran ZBETA and ZVECIN/ZVECOUT. This makes all work space
// leading dimensions independent of the mode, so that GEMMs of the same shape from different
// nodes and modes can be batched together.
//
// All index/offset arrays are host arrays, apart from lev_node_iclist, lev_node_pnonim and
// lev_node_b which, together with A and C, are device arrays. Indices stored in the flattened
// struct (IFCOL, IFROW, ICLIST) are Fortran 1-based; IOFFBETA and all *_OFFSET arrays are 0-based.
//
// The only true dependencies are between consecutive levels of the same mode (and, in the
// transposed case, between the two partner nodes scattering into the same children), so all modes
// are advanced together one level at a time. For each level a plan is built on the host: one task
// per node for the gather/scatter kernels, and the GEMMs grouped by shape. The plan for all levels
// is uploaded with a single copy, after which each level takes one gather kernel (or two scatter
// kernels in the transposed case) plus one batched GEMM per distinct shape. The GEMM groups of a
// phase are independent. ECTRANS_BUTTERFLY_NBRANCHES selects independent work (0) or stream chains (N,
// default 64 for FP32 cuBLASLt, 4 otherwise) so that they can run concurrently.
// CUDA branches share one BLAS handle: host enqueueing is serial and every GEMM
// selects its stream explicitly. This avoids allocating a workspace per handle.
// HIP retains per-stream handles because rocBLAS scratch memory is stream-specific.
//
// By default the whole sequence is captured once in a CUDA/HIP graph, together with its own copy
// of the plan, and replayed on subsequent calls, which avoids both building the plan and issuing
// the individual launches. The graph is recaptured if anything it depends on changes (A, C, the
// work space or the butterfly structure). Set ECTRANS_BUTTERFLY_GRAPHS=0 for direct execution.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

#include "hicblas.h"
#include "gemm_execution.h"
#ifdef CUDAGPU
#include "hicblas_cublaslt.cuda.h"
#endif

namespace {

// Device work space, grown as needed and kept between calls
template <typename Real> Real *get_workspace(size_t n) {
  static Real *ptr = nullptr;
  static size_t capacity = 0;
  if (n > capacity) {
    // Previous calls may still be using the work space
    HIC_CHECK(hipDeviceSynchronize());
    if (ptr) HIC_CHECK(hipFree(ptr));
    HIC_CHECK(hipMalloc(&ptr, n * sizeof(Real)));
    capacity = n;
  }
  return ptr;
}

// Pinned host staging buffer and device buffer for the plan, grown as needed and kept between
// calls. The staging buffer is only refilled once the previous upload from it has completed.
struct PlanBuffers {
  char *host = nullptr;
  char *device = nullptr;
  size_t capacity = 0;
  hipEvent_t uploaded;
  bool pending = false;
};

PlanBuffers &get_plan_buffers(size_t n) {
  static PlanBuffers buf;
  static bool initialised = false;
  if (!initialised) {
    HIC_CHECK(hipEventCreateWithFlags(&buf.uploaded, hipEventDisableTiming));
    initialised = true;
  }
  if (buf.pending) {
    HIC_CHECK(hipEventSynchronize(buf.uploaded));
    buf.pending = false;
  }
  if (n > buf.capacity) {
    // Previous calls may still be reading the device buffer
    HIC_CHECK(hipDeviceSynchronize());
    if (buf.host) HIC_CHECK(hipHostFree(buf.host));
    if (buf.device) HIC_CHECK(hipFree(buf.device));
    HIC_CHECK(hipHostMalloc((void **)&buf.host, n, 0));
    HIC_CHECK(hipMalloc(&buf.device, n));
    buf.capacity = n;
  }
  return buf;
}

void gemm(hipStream_t stream, hipblasOperation_t transa, hipblasOperation_t transb, int m,
          int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta,
          float *C, int ldc) {
#ifdef CUDAGPU
  if (legendre_cublaslt::enabled()) {
    legendre_cublaslt::gemm(stream, transa == HIPBLAS_OP_T, transb == HIPBLAS_OP_T,
                            m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    return;
  }
#endif
  HICBLAS_CHECK(hipblasSgemm(gemm_execution::blas().handle(stream), transa, transb, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc));
}

void gemm(hipStream_t stream, hipblasOperation_t transa, hipblasOperation_t transb, int m,
          int n, int k, double alpha, const double *A, int lda, const double *B, int ldb,
          double beta, double *C, int ldc) {
  HICBLAS_CHECK(hipblasDgemm(gemm_execution::blas().handle(stream), transa, transb, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc));
}

void gemm_batched(hipStream_t stream, hipblasOperation_t transa, hipblasOperation_t transb,
                  int m, int n, int k, float alpha, const float *const *A, int lda,
                  const float *const *B, int ldb, float beta, float *const *C, int ldc,
                  int batch_count) {
#ifdef CUDAGPU
  if (legendre_cublaslt::enabled()) {
    legendre_cublaslt::gemm(stream, transa == HIPBLAS_OP_T, transb == HIPBLAS_OP_T,
                            m, n, k, alpha, A, lda, B, ldb, beta,
                            const_cast<void *>(static_cast<const void *>(C)),
                            ldc, batch_count, true);
    return;
  }
#endif
  HICBLAS_CHECK(hipblasSgemmBatched(gemm_execution::blas().handle(stream), transa, transb, m, n, k, &alpha, A, lda, B, ldb, &beta,
                                    C, ldc, batch_count));
}

void gemm_batched(hipStream_t stream, hipblasOperation_t transa, hipblasOperation_t transb,
                  int m, int n, int k, double alpha, const double *const *A, int lda,
                  const double *const *B, int ldb, double beta, double *const *C, int ldc,
                  int batch_count) {
  HICBLAS_CHECK(hipblasDgemmBatched(gemm_execution::blas().handle(stream), transa, transb, m, n, k, &alpha, A, lda, B, ldb, &beta,
                                    C, ldc, batch_count));
}

void check_rank(int rank) {
  if (rank <= 0) {
    fprintf(stderr, "mult_butm: IRANK<=0 not allowed\n");
    abort();
  }
}

constexpr int block_size = 256;

// Gather/scatter task for one node. Offsets are in elements; ext_off is relative to A (gather)
// or C (scatter) if the EXTERNAL flag is set, and relative to the work space otherwise.
enum TaskFlags { EXTERNAL = 1, ACCUMULATE = 2 };

struct NodeTask {
  int64_t ext_off;  // Gather source / scatter destination (row 0, field 0)
  int64_t beta_off; // This node's beta slot
  int64_t vec_off;  // This node's vec region
  int icols, rank, iclist_off;
  int off_l, split, off_r, limit; // Scatter only: see scatter_kernel
  int flags;
};

// Small node counts otherwise leave most SMs idle while a single block copies a
// large node. Spread each node over several blocks, targeting one full wave of
// threads across the device. Each block owns a disjoint subset of field elements.
dim3 node_task_grid(const std::vector<NodeTask> &tasks, int n_flds) {
  int tiles = 1;
#ifdef CUDAGPU
  static const int target_blocks = [] {
    int device, sms, threads;
    HIC_CHECK(cudaGetDevice(&device));
    HIC_CHECK(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device));
    HIC_CHECK(cudaDeviceGetAttribute(&threads, cudaDevAttrMaxThreadsPerMultiProcessor, device));
    return sms * std::max(1, threads / block_size);
  }();
  int64_t max_elements = 0;
  for (const auto &task : tasks)
    max_elements = std::max(max_elements, int64_t(task.icols) * n_flds);
  const int max_tiles = int(std::min<int64_t>(64, (max_elements + block_size - 1) / block_size));
  const int count = int(tasks.size());
  tiles = std::max(1, std::min(max_tiles, (target_blocks + count - 1) / count));
#endif
  return dim3(unsigned(tasks.size()), unsigned(tiles));
}

// Gather the columns of each node, permuted by ICLIST, from the source. The first rank columns go
// to beta (the "skeleton" columns), the rest go to vec (to be multiplied by PNONIM). One block per
// node, split across blocks when there are few nodes. The source is A (row stride lda) for EXTERNAL tasks, and the work space (row stride n_flds)
// otherwise.
template <typename Real>
__global__ void gather_kernel(const NodeTask *tasks, int n_flds, const int *iclist_all,
                              const Real *A, int lda, Real *work) {
  const NodeTask t = tasks[blockIdx.x];
  const int *iclist = iclist_all + t.iclist_off;
  bool external = t.flags & EXTERNAL;
  const Real *src = external ? A + t.ext_off : work + t.ext_off;
  int src_row_stride = external ? lda : n_flds;

  for (int e = threadIdx.x + blockIdx.y * blockDim.x; e < t.icols * n_flds;
       e += blockDim.x * gridDim.y) {
    int f = e % n_flds;
    int jn = e / n_flds;
    Real v = src[(size_t)(iclist[jn] - 1) * src_row_stride + f];
    Real *dst = work + (jn < t.rank ? t.beta_off : t.vec_off);
    dst[(size_t)jn * n_flds + f] = v;
  }
}

// Scatter the columns of each node, permuted by ICLIST, into the destination (transposed case).
// Each node is split across disjoint blocks as in gather_kernel. Column jn is taken from beta if jn < rank and from vec otherwise.
// Destination row i = ICLIST(jn)-1 maps to row off_l + i if i < split, to row off_r + (i - split)
// if split <= i < limit, and is dropped otherwise. With ACCUMULATE the values are added to the
// destination, otherwise they overwrite it. The destination is C (row stride ldc) for EXTERNAL
// tasks, and the work space (row stride n_flds) otherwise.
template <typename Real>
__global__ void scatter_kernel(const NodeTask *tasks, int n_flds, const int *iclist_all, Real *C,
                               int ldc, Real *work) {
  const NodeTask t = tasks[blockIdx.x];
  const int *iclist = iclist_all + t.iclist_off;
  bool external = t.flags & EXTERNAL;
  bool accumulate = t.flags & ACCUMULATE;
  Real *dst = external ? C + t.ext_off : work + t.ext_off;
  int dst_row_stride = external ? ldc : n_flds;

  for (int e = threadIdx.x + blockIdx.y * blockDim.x; e < t.icols * n_flds;
       e += blockDim.x * gridDim.y) {
    int f = e % n_flds;
    int jn = e / n_flds;
    int i = iclist[jn] - 1;
    if (i >= t.limit) continue;
    int row = i < t.split ? t.off_l + i : t.off_r + (i - t.split);

    const Real *src = work + (jn < t.rank ? t.beta_off : t.vec_off);
    Real v = src[(size_t)jn * n_flds + f];
    Real *p = dst + (size_t)row * dst_row_stride + f;
    *p = accumulate ? *p + v : v;
  }
}

// A group of GEMMs with identical shape, leading dimensions and operations, executed as one
// batched GEMM. Pointer arrays are filled on the host and uploaded with the rest of the plan;
// ptr_off is the position of the group's A array in the uploaded pointer list, followed by its B
// and C arrays.
template <typename Real> struct GemmGroup {
  hipblasOperation_t transa, transb;
  int m, n, k, lda, ldb, ldc;
  Real beta;
  std::vector<const Real *> a, b;
  std::vector<Real *> c;
  size_t ptr_off = 0;
};

// The GEMMs of one phase of a level, grouped by (m, n, k). Within a phase all GEMMs have the same
// operations and the leading dimensions follow from the shape.
template <typename Real> struct GemmPhase {
  std::map<std::tuple<int, int, int>, size_t> index;
  std::vector<GemmGroup<Real>> groups;

  void add(hipblasOperation_t transa, hipblasOperation_t transb, int m, int n, int k,
           const Real *a, int lda, const Real *b, int ldb, Real beta, Real *c, int ldc) {
    auto key = std::make_tuple(m, n, k);
    auto it = index.find(key);
    if (it == index.end()) {
      it = index.emplace(key, groups.size()).first;
      GemmGroup<Real> g;
      g.transa = transa;
      g.transb = transb;
      g.m = m;
      g.n = n;
      g.k = k;
      g.lda = lda;
      g.ldb = ldb;
      g.ldc = ldc;
      g.beta = beta;
      groups.push_back(g);
    }
    GemmGroup<Real> &g = groups[it->second];
    g.a.push_back(a);
    g.b.push_back(b);
    g.c.push_back(c);
  }

  void launch_group(int i, hipStream_t stream, void *const *d_ptrs) const {
    const GemmGroup<Real> &g = groups[i];
    const int count = int(g.c.size());
    if (count == 1) {
      gemm(stream, g.transa, g.transb, g.m, g.n, g.k, Real(1),
           g.a[0], g.lda, g.b[0], g.ldb, g.beta, g.c[0], g.ldc);
    } else {
      void *const *p = d_ptrs + g.ptr_off;
      gemm_batched(stream, g.transa, g.transb, g.m, g.n, g.k, Real(1),
                   reinterpret_cast<const Real *const *>(p), g.lda,
                   reinterpret_cast<const Real *const *>(p + count), g.ldb, g.beta,
                   reinterpret_cast<Real *const *>(p + 2 * count), g.ldc, count);
    }
  }

  void launch(gemm_execution::Scheduler &scheduler, hipStream_t stream,
              void *const *d_ptrs) const {
    scheduler.run(int(groups.size()), stream,
                  [&](int i, hipStream_t branch) { launch_group(i, branch, d_ptrs); });
  }
};

// Plan for one level (across all modes). In the normal case the gather tasks are all in
// tasks[0], followed by gemms[0] (ZBETA += PNONIM * ZVECIN) and gemms[1] (PVECOUT = B * ZBETA, for
// modes whose last level this is). In the transposed case gemms[0] (ZBETA = B^T * PVECIN, for
// modes whose first level this is) and gemms[1] (ZVECOUT = PNONIM^T * ZBETA) come first, followed
// by the scatter tasks that assign (tasks[0]) and then those that accumulate (tasks[1]).
template <typename Real> struct LevelPlan {
  std::vector<NodeTask> tasks[2];
  size_t task_off[2] = {0, 0};
  GemmPhase<Real> gemms[2];
};

#ifdef CUDAGPU
template <typename Real>
void prepare_cublaslt(const std::vector<LevelPlan<Real>> &plan) {
  if constexpr (std::is_same<Real, float>::value) {
    if (!legendre_cublaslt::enabled()) return;
    for (const auto &level : plan)
      for (int phase = 0; phase < 2; ++phase)
        for (const auto &group : level.gemms[phase].groups) {
          const int count = int(group.c.size());
          if (count == 0) continue;
          legendre_cublaslt::prepare(
              group.transa == HIPBLAS_OP_T, group.transb == HIPBLAS_OP_T,
              group.m, group.n, group.k,
              group.lda, group.ldb, group.ldc, count, count > 1,
              group.a[0], group.b[0], group.c[0]);
        }
  }
}
#endif

// Graph cache. There is at most one graph per key; everything else the graph depends on is kept
// in GraphDeps, and if any of it changes the graph is recaptured.
struct GraphKey {
  char transpose;
  int n_flds, lda, ldc, n_modes;
  const void *iclist; // Identifies the butterfly structure

  bool operator<(const GraphKey &o) const {
    return std::tie(transpose, n_flds, lda, ldc, n_modes, iclist) <
           std::tie(o.transpose, o.n_flds, o.lda, o.ldc, o.n_modes, o.iclist);
  }
};

struct GraphDeps {
  const void *A, *C, *work, *pnonim, *b;
  uint64_t hash; // Of the per-mode sizes and offsets, and the host node array addresses

  bool operator==(const GraphDeps &o) const {
    return std::tie(A, C, work, pnonim, b, hash) ==
           std::tie(o.A, o.C, o.work, o.pnonim, o.b, o.hash);
  }
};

struct GraphEntry {
  GraphDeps deps;
  size_t work_len;
  char *d_plan; // The graph's own copy of the plan
  hipGraphExec_t exec;
};

template <typename Real> std::map<GraphKey, GraphEntry> &get_graph_cache() {
  static std::map<GraphKey, GraphEntry> cache;
  return cache;
}

// FNV-1a
struct Hasher {
  uint64_t h = 14695981039346656037ull;
  void add(const void *p, size_t n) {
    const unsigned char *c = static_cast<const unsigned char *>(p);
    for (size_t i = 0; i < n; ++i) {
      h ^= c[i];
      h *= 1099511628211ull;
    }
  }
  template <typename T> void add_value(T v) { add(&v, sizeof(v)); }
};

} // namespace

template <typename Real>
void mult_butm(char transpose, int start_mode, int end_mode, int n_flds, const int *order,
  const int *levels, const int *betalen_max, const int *lev_offset, const int *lev_ij,
  const int *lev_ik, const int *lev_ibetalen, const int *lev_node_offset, const int *lev_node_ifcol,
  const int *lev_node_ilcol, const int *lev_node_ifrow, const int *lev_node_ilrow,
  const int *lev_node_icols, const int *lev_node_irows, const int *lev_node_irank,
  const int *lev_node_ioffbeta, const int *lev_node_iclist_offset, const int *lev_node_iclist,
  const int *lev_node_pnonim_offset, const Real *lev_node_pnonim, const int *lev_node_b_offset,
  const Real *lev_node_b, const Real *A, int lda, const int64_t *a_offsets, Real *C, int ldc,
  const int64_t *c_offsets, hipStream_t stream) {

  bool is_transposed = (transpose == 'T' || transpose == 't');

  int n_modes = end_mode - start_mode + 1;
  if (n_modes <= 0 || n_flds <= 0) return;

  int default_branches = 4;
#ifdef CUDAGPU
  if constexpr (std::is_same<Real, float>::value)
    if (legendre_cublaslt::enabled()) default_branches = 64;
#endif
  static gemm_execution::Scheduler scheduler(gemm_execution::options(true, default_branches));

  const hipblasOperation_t N = HIPBLAS_OP_N, T = HIPBLAS_OP_T;
  const int nf = n_flds;

  // Work space for each mode: two "beta" buffers (ZBETA(:,:,0:1)) followed by the vec regions
  // (ZVECIN in the normal case, ZVECOUT in the transposed case) of the nodes of one level. Each
  // node gets its own ICOLS x n_flds vec region, so the vec part is sized for the level with the
  // largest total ICOLS.
  std::vector<size_t> work_offsets(n_modes + 1, 0);
  int max_levels = 0;
  auto compute_work_offsets = [&]() {
    for (int m = start_mode; m <= end_mode; ++m) {
      size_t vec_len = 0;
      for (int lev = 0; lev <= levels[m]; ++lev) {
        int l = lev_offset[m] + lev;
        size_t lev_len = 0;
        for (int inode = 0; inode < lev_ij[l] * lev_ik[l]; ++inode) {
          lev_len += lev_node_icols[lev_node_offset[l] + inode];
        }
        vec_len = std::max(vec_len, lev_len);
      }
      work_offsets[m + 1 - start_mode] = work_offsets[m - start_mode] + (2 * (size_t)betalen_max[m] + vec_len) * nf;
      max_levels = std::max(max_levels, levels[m]);
    }
  };

  // Build the plan for every level, for the work space at work
  auto build_plan = [&](Real *work) {
    std::vector<LevelPlan<Real>> plan(max_levels + 1);
    for (int m = start_mode; m <= end_mode; ++m) {
      int nlevels = levels[m];
      int lbeta = betalen_max[m];
      size_t beta_buf[2] = {work_offsets[m - start_mode], work_offsets[m - start_mode] + (size_t)(lbeta) * nf};
      size_t vec_base = work_offsets[m - start_mode] + (size_t)2 * lbeta * nf;

      // Flat index of NODE(j,k) on level lev of this mode (j and k are 1-based, as in Fortran)
      auto node_index = [&](int lev, int j, int k) {
        int l = lev_offset[m] + lev;
        return lev_node_offset[l] + (k - 1) * lev_ij[l] + (j - 1);
      };

      for (int lev = 0; lev <= nlevels; ++lev) {
        int l = lev_offset[m] + lev;
        LevelPlan<Real> &lp = plan[lev];
        size_t beta = beta_buf[lev % 2];
        size_t beta_prev = beta_buf[(lev + 1) % 2]; // Buffer for level lev-1
        size_t vec_off = vec_base;
        for (int j = 1; j <= lev_ij[l]; ++j) {
          for (int k = 1; k <= lev_ik[l]; ++k) {
            int idx = node_index(lev, j, k);
            int rank = lev_node_irank[idx];
            int icols = lev_node_icols[idx];
            int irows = lev_node_irows[idx];
            int n = icols - rank;
            check_rank(rank);

            NodeTask t = {};
            t.icols = icols;
            t.rank = rank;
            t.iclist_off = lev_node_iclist_offset[idx];
            t.beta_off = beta + (size_t)lev_node_ioffbeta[idx] * nf;
            t.vec_off = vec_off;
            vec_off += (size_t)icols * nf;

            Real *beta_ptr = work + t.beta_off;
            Real *vec_ns_ptr = work + t.vec_off + (size_t)rank * nf; // Non-skeleton columns of vec
            const Real *pnonim = lev_node_pnonim + lev_node_pnonim_offset[idx];
            const Real *b = lev_node_b + lev_node_b_offset[idx];

            if (is_transposed) {
              if (lev == nlevels) {
                // ZBETA = B^T * PVECIN(IFR:ILR,:), i.e. ZBETA^T = PVECIN(IFR:ILR,:)^T * B
                lp.gemms[0].add(N, N, nf, rank, irows,
                                A + a_offsets[m] + (size_t)(lev_node_ifrow[idx] - 1) * lda, lda,
                                b, irows, (Real)0.0, beta_ptr, nf);
              }

              if (n > 0) {
                // ZVECOUT(IRANK+1:ICOLS,:) = PNONIM^T * ZBETA, i.e. ZVECOUT^T = ZBETA^T * PNONIM
                lp.gemms[1].add(N, N, nf, n, rank, beta_ptr, nf, pnonim, rank, (Real)0.0,
                                vec_ns_ptr, nf);
              }

              if (lev == 0) {
                // Scatter into PVECOUT(IFR:ILR,:)
                t.flags = EXTERNAL;
                t.ext_off = c_offsets[m] + (int64_t)(lev_node_ifcol[idx] - 1) * ldc;
                t.off_l = 0;
                t.split = icols;
                t.off_r = 0;
                t.limit = icols;
                lp.tasks[0].push_back(t);
              } else {
                // Scatter into the beta slots of the two children on level lev-1. Odd j assigns,
                // even j (the partner sharing the same children) accumulates afterwards.
                int jc = (j + 1) / 2;
                int il = node_index(lev - 1, jc, 2 * k - 1);
                int irankl = lev_node_irank[il];
                int irankr = 0, btstr = 0;
                if (2 * k <= lev_ik[l - 1]) {
                  int ir = node_index(lev - 1, jc, 2 * k);
                  irankr = lev_node_irank[ir];
                  btstr = lev_node_ioffbeta[ir];
                }
                t.ext_off = beta_prev;
                t.off_l = lev_node_ioffbeta[il];
                t.split = irankl;
                t.off_r = btstr;
                t.limit = irankl + irankr;
                if (j % 2 == 0) {
                  t.flags = ACCUMULATE;
                  lp.tasks[1].push_back(t);
                } else {
                  lp.tasks[0].push_back(t);
                }
              }
            } else {
              if (lev == 0) {
                // Gather from PVECIN(IFR:,:)
                t.flags = EXTERNAL;
                t.ext_off = a_offsets[m] + (int64_t)(lev_node_ifcol[idx] - 1) * lda;
              } else {
                // Gather from the beta slots of the children on level lev-1, which are contiguous
                // starting from the left child
                int il = node_index(lev - 1, (j + 1) / 2, 2 * k - 1);
                t.ext_off = beta_prev + (size_t)lev_node_ioffbeta[il] * nf;
              }
              lp.tasks[0].push_back(t);

              if (n > 0) {
                // ZBETA += PNONIM * ZVECIN(IRANK+1:ICOLS,:), i.e.
                // ZBETA^T += ZVECIN(IRANK+1:ICOLS,:)^T * PNONIM^T
                lp.gemms[0].add(N, T, nf, rank, n, vec_ns_ptr, nf, pnonim, rank, (Real)1.0,
                                beta_ptr, nf);
              }

              if (lev == nlevels) {
                // PVECOUT(IFR:ILR,:) = B * ZBETA, i.e. PVECOUT(IFR:ILR,:)^T = ZBETA^T * B^T
                lp.gemms[1].add(N, T, nf, irows, rank, beta_ptr, nf, b, irows, (Real)0.0,
                                C + c_offsets[m] + (size_t)(lev_node_ifrow[idx] - 1) * ldc, ldc);
              }
            }
          }
        }
      }
    }
    return plan;
  };

  // Pack the tasks and GEMM pointer arrays of all levels into one buffer, to be uploaded in one
  // go. Returns the size of the tasks part, which is followed by the pointers.
  auto pack_plan = [&](std::vector<LevelPlan<Real>> &plan, std::vector<char> &host) {
    size_t n_tasks = 0, n_ptrs = 0;
    for (LevelPlan<Real> &lp : plan) {
      for (int i = 0; i < 2; ++i) {
        lp.task_off[i] = n_tasks;
        n_tasks += lp.tasks[i].size();
        for (GemmGroup<Real> &g : lp.gemms[i].groups) {
          g.ptr_off = n_ptrs;
          n_ptrs += 3 * g.c.size();
        }
      }
    }
    size_t task_bytes = n_tasks * sizeof(NodeTask);
    host.resize(task_bytes + n_ptrs * sizeof(void *));
    NodeTask *h_tasks = reinterpret_cast<NodeTask *>(host.data());
    const void **h_ptrs = reinterpret_cast<const void **>(host.data() + task_bytes);
    for (const LevelPlan<Real> &lp : plan) {
      for (int i = 0; i < 2; ++i) {
        std::copy(lp.tasks[i].begin(), lp.tasks[i].end(), h_tasks + lp.task_off[i]);
        for (const GemmGroup<Real> &g : lp.gemms[i].groups) {
          size_t count = g.c.size();
          std::copy(g.a.begin(), g.a.end(), h_ptrs + g.ptr_off);
          std::copy(g.b.begin(), g.b.end(), h_ptrs + g.ptr_off + count);
          std::copy(g.c.begin(), g.c.end(), h_ptrs + g.ptr_off + 2 * count);
        }
      }
    }
    return task_bytes;
  };

  // Keep the numerical dependencies identical for direct execution and child graphs.
  auto visit_plan = [&](const std::vector<LevelPlan<Real>> &plan, auto tasks, auto phase) {
    if (is_transposed) {
      for (int lev = max_levels; lev >= 0; --lev) {
        const auto &lp = plan[lev];
        phase(lp.gemms[0]);
        phase(lp.gemms[1]);
        tasks(lp, 0, true);
        tasks(lp, 1, true);
      }
    } else {
      for (const auto &lp : plan) {
        tasks(lp, 0, false);
        phase(lp.gemms[0]);
        phase(lp.gemms[1]);
      }
    }
  };

  auto launch_tasks = [&](const LevelPlan<Real> &lp, int i, bool scatter,
                          const NodeTask *d_tasks, Real *work, hipStream_t stream) {
    if (lp.tasks[i].empty()) return;
    const dim3 grid = node_task_grid(lp.tasks[i], nf);
    if (scatter)
      scatter_kernel<<<grid, block_size, 0, stream>>>(d_tasks + lp.task_off[i], nf,
                                                     lev_node_iclist, C, ldc, work);
    else
      gather_kernel<<<grid, block_size, 0, stream>>>(d_tasks + lp.task_off[i], nf,
                                                    lev_node_iclist, A, lda, work);
    HIC_CHECK(hipGetLastError());
  };

  auto prepare_execution = [&](const std::vector<LevelPlan<Real>> &plan) {
    int max_groups = 0;
    for (const auto &level : plan)
      for (const auto &phase : level.gemms)
        max_groups = std::max(max_groups, int(phase.groups.size()));
    scheduler.prepare(max_groups, stream, [&](hipStream_t branch) {
#ifdef CUDAGPU
      if constexpr (std::is_same<Real, float>::value)
        if (legendre_cublaslt::enabled()) return;
#endif
      gemm_execution::blas().prepare<Real>(branch, scheduler.config.graphs && scheduler.config.branches == 0);
    });
#ifdef CUDAGPU
    if (scheduler.config.graphs) prepare_cublaslt(plan);
#endif
  };

  auto execute_plan = [&](const std::vector<LevelPlan<Real>> &plan, char *d_plan,
                          size_t task_bytes, Real *work, hipStream_t stream) {
    const NodeTask *d_tasks = reinterpret_cast<const NodeTask *>(d_plan);
    void *const *d_ptrs = reinterpret_cast<void *const *>(d_plan + task_bytes);
    visit_plan(plan,
        [&](const auto &lp, int i, bool scatter) { launch_tasks(lp, i, scatter, d_tasks, work, stream); },
        [&](const auto &phase) { phase.launch(scheduler, stream, d_ptrs); });
  };

  const bool use_graphs = scheduler.config.graphs;

  if (!use_graphs) {
    // Build, upload and execute the plan directly on the caller's stream
    compute_work_offsets();
    Real *work = get_workspace<Real>(work_offsets[n_modes]);
    std::vector<LevelPlan<Real>> plan = build_plan(work);
    prepare_execution(plan);
    std::vector<char> host;
    size_t task_bytes = pack_plan(plan, host);
    PlanBuffers &buf = get_plan_buffers(host.size());
    if (!host.empty()) {
      std::copy(host.begin(), host.end(), buf.host);
      HIC_CHECK(hipMemcpyAsync(buf.device, buf.host, host.size(), hipMemcpyHostToDevice, stream));
      HIC_CHECK(hipEventRecord(buf.uploaded, stream));
      buf.pending = true;
    }
    execute_plan(plan, buf.device, task_bytes, work, stream);
    return;
  }

  // Everything, apart from the key, that a captured graph depends on. The sizes of the butterfly
  // structure are hashed per mode and per level only, since hashing all the node arrays on every
  // call would be relatively expensive; the device arrays identify the structure itself.
  auto make_deps = [&](const Real *work) {
    int n_lev_total = lev_offset[end_mode] + levels[end_mode] + 1 - lev_offset[start_mode];
    Hasher h;
    h.add(order, n_modes * sizeof(int));
    h.add(levels, n_modes * sizeof(int));
    h.add(betalen_max, n_modes * sizeof(int));
    h.add(lev_offset, n_modes * sizeof(int));
    h.add(a_offsets, n_modes * sizeof(int64_t));
    h.add(c_offsets, n_modes * sizeof(int64_t));
    h.add(lev_ij, n_lev_total * sizeof(int));
    h.add(lev_ik, n_lev_total * sizeof(int));
    h.add(lev_node_offset, n_lev_total * sizeof(int));
    return GraphDeps{A, C, work, lev_node_pnonim, lev_node_b, h.h};
  };

  GraphKey key = {is_transposed ? 'T' : 'N', n_flds, lda, ldc, n_modes, lev_node_iclist};
  auto &cache = get_graph_cache<Real>();
  auto it = cache.find(key);
  if (it != cache.end()) {
    GraphEntry &entry = it->second;
    Real *work = get_workspace<Real>(entry.work_len);
    if (make_deps(work) == entry.deps) {
      HIC_CHECK(hipGraphLaunch(entry.exec, stream));
      return;
    }
    fprintf(stderr, "WARNING mult_butm: graph dependencies changed, recapturing (slow)\n");
    HIC_CHECK(hipDeviceSynchronize());
    HIC_CHECK(hipGraphExecDestroy(entry.exec));
    if (entry.d_plan) HIC_CHECK(hipFree(entry.d_plan));
    cache.erase(it);
  }

  // Build the plan and give the graph its own copy of it on the device
  compute_work_offsets();
  GraphEntry entry;
  entry.work_len = work_offsets[n_modes];
  Real *work = get_workspace<Real>(entry.work_len);
  std::vector<LevelPlan<Real>> plan = build_plan(work);
  prepare_execution(plan);
  std::vector<char> host;
  size_t task_bytes = pack_plan(plan, host);
  entry.d_plan = nullptr;
  if (!host.empty()) {
    HIC_CHECK(hipMalloc(&entry.d_plan, host.size()));
    HIC_CHECK(hipMemcpy(entry.d_plan, host.data(), host.size(), hipMemcpyHostToDevice));
  }
  entry.deps = make_deps(work);

  hipGraph_t graph;
  if (scheduler.config.branches == 0) {
    HIC_CHECK(hipGraphCreate(&graph, 0));
    std::vector<hipGraphNode_t> parents;
    const NodeTask *d_tasks = reinterpret_cast<const NodeTask *>(entry.d_plan);
    void *const *d_ptrs = reinterpret_cast<void *const *>(entry.d_plan + task_bytes);
    visit_plan(plan, [&](const auto &lp, int i, bool scatter) {
      if (lp.tasks[i].empty()) return;
      auto node = scheduler.append_serial(graph, parents, scheduler.capture_stream(),
          [&](hipStream_t branch) { launch_tasks(lp, i, scatter, d_tasks, work, branch); });
      parents = {node};
    }, [&](const auto &phase) {
      scheduler.append(graph, parents, int(phase.groups.size()),
          [&](int i, hipStream_t branch) { phase.launch_group(i, branch, d_ptrs); });
    });
  } else {
    // Retain the single captured fork/join graph for bounded stream chains.
    HIC_CHECK(hipStreamBeginCapture(scheduler.capture_stream(), hipStreamCaptureModeGlobal));
    execute_plan(plan, entry.d_plan, task_bytes, work, scheduler.capture_stream());
    HIC_CHECK(hipStreamEndCapture(scheduler.capture_stream(), &graph));
  }
  HIC_CHECK(hipGraphInstantiate(&entry.exec, graph, NULL, NULL, 0));
  HIC_CHECK(hipGraphDestroy(graph));
  cache.emplace(key, entry);

  HIC_CHECK(hipGraphLaunch(entry.exec, stream));
}

// -------------------------------------------------------------------------------------------------
// Fortran-callable wrappers for the mult_butm function
// -------------------------------------------------------------------------------------------------

extern "C" {
void mult_butm_sp(
  char transpose, int start_mode, int end_mode, int n_flds, const int *order, const int *levels,
  const int *betalen_max, const int *lev_offset, const int *lev_ij, const int *lev_ik,
  const int *lev_ibetalen, const int *lev_node_offset, const int *lev_node_ifcol,
  const int *lev_node_ilcol, const int *lev_node_ifrow, const int *lev_node_ilrow,
  const int *lev_node_icols, const int *lev_node_irows, const int *lev_node_irank,
  const int *lev_node_ioffbeta, const int *lev_node_iclist_offset, const int *lev_node_iclist,
  const int *lev_node_pnonim_offset, const float *lev_node_pnonim, const int *lev_node_b_offset,
  const float *lev_node_b, const float *A, int lda, const int64_t *a_offsets, float *C, int ldc,
  const int64_t *c_offsets, const long *stream) {

  mult_butm<float>(transpose, start_mode, end_mode, n_flds, order, levels, betalen_max,
    lev_offset, lev_ij, lev_ik, lev_ibetalen,
    lev_node_offset, lev_node_ifcol, lev_node_ilcol,
    lev_node_ifrow, lev_node_ilrow, lev_node_icols,
    lev_node_irows, lev_node_irank, lev_node_ioffbeta,
    lev_node_iclist_offset, lev_node_iclist, lev_node_pnonim_offset,
    lev_node_pnonim, lev_node_b_offset, lev_node_b,
    A, lda, a_offsets, C, ldc, c_offsets, reinterpret_cast<hipStream_t>(*stream));
}

void mult_butm_dp(
  char transpose, int start_mode, int end_mode, int n_flds, const int *order, const int *levels,
  const int *betalen_max, const int *lev_offset, const int *lev_ij, const int *lev_ik,
  const int *lev_ibetalen, const int *lev_node_offset, const int *lev_node_ifcol,
  const int *lev_node_ilcol, const int *lev_node_ifrow, const int *lev_node_ilrow,
  const int *lev_node_icols, const int *lev_node_irows, const int *lev_node_irank,
  const int *lev_node_ioffbeta, const int *lev_node_iclist_offset, const int *lev_node_iclist,
  const int *lev_node_pnonim_offset, const double *lev_node_pnonim, const int *lev_node_b_offset,
  const double *lev_node_b, const double *A, int lda, const int64_t *a_offsets, double *C, int ldc,
  const int64_t *c_offsets, const long *stream) {

  mult_butm<double>(transpose, start_mode, end_mode, n_flds, order, levels, betalen_max,
    lev_offset, lev_ij, lev_ik, lev_ibetalen,
    lev_node_offset, lev_node_ifcol, lev_node_ilcol,
    lev_node_ifrow, lev_node_ilrow, lev_node_icols,
    lev_node_irows, lev_node_irank, lev_node_ioffbeta,
    lev_node_iclist_offset, lev_node_iclist, lev_node_pnonim_offset,
    lev_node_pnonim, lev_node_b_offset, lev_node_b,
    A, lda, a_offsets, C, ldc, c_offsets, reinterpret_cast<hipStream_t>(*stream));
}
}
