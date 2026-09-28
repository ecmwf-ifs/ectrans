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
// Unlike the CPU version, the input (A) and output (C) arrays are stored with the field index
// fastest, as for the other GPU Legendre GEMMs, i.e. element (row, field) of mode m is at
// A[a_offsets[m] + field + row * lda]. The internal work arrays (beta and vec) keep the CPU
// layout (row fastest, leading dimensions IBETALEN_MAX and ICOLS respectively) so the logic
// below maps directly onto the Fortran.
//
// All index/offset arrays are host arrays, apart from lev_node_iclist, lev_node_pnonim and
// lev_node_b which, together with A and C, are device arrays. Indices stored in the flattened
// struct (IFCOL, IFROW, ICLIST) are Fortran 1-based; IOFFBETA and all *_OFFSET arrays are 0-based.
//
// Modes are distributed round-robin over a pool of worker streams (MULT_BUTM_NSTREAMS, default
// 4), which are forked from and joined back onto the caller's stream. Each stream has its own
// work space, and each node on a level its own vec region within it, so that the only true
// dependencies within a mode are between consecutive levels (and, in the transposed case,
// between the two partner nodes scattering into the same children).

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "hicblas.h"

namespace {

// Pool of worker streams, each with its own BLAS handle, plus the events used to fork from and
// join onto the caller's stream
struct StreamPool {
  int n = 0;
  std::vector<hipStream_t> streams;
  std::vector<hipblasHandle_t> handles;
  std::vector<hipEvent_t> join_events;
  hipEvent_t fork_event;
};

StreamPool &get_stream_pool() {
  static StreamPool pool;
  if (pool.n == 0) {
    int n = 4;
    if (const char *env = std::getenv("MULT_BUTM_NSTREAMS")) n = std::max(1, std::atoi(env));
    pool.streams.resize(n);
    pool.handles.resize(n);
    pool.join_events.resize(n);
    for (int i = 0; i < n; ++i) {
      HIC_CHECK(hipStreamCreateWithFlags(&pool.streams[i], hipStreamNonBlocking));
      HICBLAS_CHECK(hipblasCreate(&pool.handles[i]));
      HICBLAS_CHECK(hipblasSetStream(pool.handles[i], pool.streams[i]));
      HIC_CHECK(hipEventCreateWithFlags(&pool.join_events[i], hipEventDisableTiming));
    }
    HIC_CHECK(hipEventCreateWithFlags(&pool.fork_event, hipEventDisableTiming));
    pool.n = n;
  }
  return pool;
}

// Device work space, grown as needed and kept between calls
template <typename Real> Real *get_workspace(size_t n) {
  static Real *ptr = nullptr;
  static size_t capacity = 0;
  if (n > capacity) {
    if (ptr) HIC_CHECK(hipFree(ptr));
    HIC_CHECK(hipMalloc(&ptr, n * sizeof(Real)));
    capacity = n;
  }
  return ptr;
}

void gemm(hipblasHandle_t handle, hipblasOperation_t transa, hipblasOperation_t transb, int m,
          int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta,
          float *C, int ldc) {
  HICBLAS_CHECK(hipblasSgemm(handle, transa, transb, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc));
}

void gemm(hipblasHandle_t handle, hipblasOperation_t transa, hipblasOperation_t transb, int m,
          int n, int k, double alpha, const double *A, int lda, const double *B, int ldb,
          double beta, double *C, int ldc) {
  HICBLAS_CHECK(hipblasDgemm(handle, transa, transb, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc));
}

void check_rank(int rank) {
  if (rank <= 0) {
    fprintf(stderr, "mult_butm: IRANK<=0 not allowed\n");
    abort();
  }
}

constexpr int block_size = 256;

inline int num_blocks(int n) { return (n + block_size - 1) / block_size; }

// Gather the columns of a node, permuted by ICLIST, from src. The first rank columns go to beta
// (the "skeleton" columns), the rest go to vec (to be multiplied by PNONIM). Element (row, field)
// of src is at src[row * src_row_stride + field * src_fld_stride].
template <typename Real>
__global__ void gather_kernel(int icols, int rank, int n_flds, const int *iclist,
                              const Real *src, int src_row_stride, int src_fld_stride,
                              Real *beta, int lbeta, Real *vec, int ld_vec) {
  int t = blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= icols * n_flds) return;
  int f = t % n_flds;
  int jn = t / n_flds;

  Real v = src[(iclist[jn] - 1) * src_row_stride + f * src_fld_stride];
  if (jn < rank) {
    beta[jn + f * lbeta] = v;
  } else {
    vec[jn + f * ld_vec] = v;
  }
}

// Scatter the columns of a node, permuted by ICLIST, into dst (transposed case). Column jn is
// taken from beta if jn < rank and from vec otherwise. Destination row i = ICLIST(jn)-1 maps to
// row off_l + i if i < split, to row off_r + (i - split) if split <= i < limit, and is dropped
// otherwise. With accumulate the values are added to dst, otherwise they overwrite it. Element
// (row, field) of dst is at dst[row * dst_row_stride + field * dst_fld_stride].
template <typename Real>
__global__ void scatter_kernel(int icols, int rank, int n_flds, const int *iclist,
                               const Real *beta, int lbeta, const Real *vec, int ld_vec,
                               Real *dst, int dst_row_stride, int dst_fld_stride,
                               int off_l, int split, int off_r, int limit, bool accumulate) {
  int t = blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= icols * n_flds) return;
  int f = t % n_flds;
  int jn = t / n_flds;

  int i = iclist[jn] - 1;
  if (i >= limit) return;
  int row = i < split ? off_l + i : off_r + (i - split);

  Real v = jn < rank ? beta[jn + f * lbeta] : vec[jn + f * ld_vec];
  Real *p = dst + row * dst_row_stride + f * dst_fld_stride;
  *p = accumulate ? *p + v : v;
}

} // namespace

template <typename Real>
void mult_butm(char transpose, int n_modes, int n_flds, const int *order, const int *levels,
  const int *betalen_max, const int *lev_offset, const int *lev_ij, const int *lev_ik,
  const int *lev_ibetalen, const int *lev_node_offset, const int *lev_node_ifcol,
  const int *lev_node_ilcol, const int *lev_node_ifrow, const int *lev_node_ilrow,
  const int *lev_node_icols, const int *lev_node_irows, const int *lev_node_irank,
  const int *lev_node_ioffbeta, const int *lev_node_iclist_offset, const int *lev_node_iclist,
  const int *lev_node_pnonim_offset, const Real *lev_node_pnonim, const int *lev_node_b_offset,
  const Real *lev_node_b, const Real *A, int lda, const int64_t *a_offsets, Real *C, int ldc,
  const int64_t *c_offsets, hipStream_t stream) {

  bool is_transposed = (transpose == 'T' || transpose == 't');

  if (n_modes <= 0 || n_flds <= 0) return;

  StreamPool &pool = get_stream_pool();

  // Work space for each worker stream: two "beta" buffers (ZBETA(:,:,0:1)) followed by the vec
  // regions (ZVECIN in the normal case, ZVECOUT in the transposed case) of the nodes of one
  // level. Each node gets its own ICOLS x n_flds vec region, so the vec part is sized for the
  // level with the largest total ICOLS. Modes on the same stream run one after the other and
  // share that stream's work space, which is sized for the largest mode.
  size_t lbeta_max = 0, vec_len_max = 0;
  for (int m = 0; m < n_modes; ++m) {
    size_t vec_len = 0;
    for (int lev = 0; lev <= levels[m]; ++lev) {
      int l = lev_offset[m] + lev;
      size_t lev_len = 0;
      for (int inode = 0; inode < lev_ij[l] * lev_ik[l]; ++inode) {
        lev_len += lev_node_icols[lev_node_offset[l] + inode];
      }
      vec_len = std::max(vec_len, lev_len);
    }
    lbeta_max = std::max(lbeta_max, (size_t)betalen_max[m]);
    vec_len_max = std::max(vec_len_max, vec_len);
  }
  size_t work_len = (2 * lbeta_max + vec_len_max) * n_flds;
  Real *work = get_workspace<Real>(work_len * pool.n);

  // Fork the worker streams from the caller's stream
  HIC_CHECK(hipEventRecord(pool.fork_event, stream));
  for (int i = 0; i < pool.n; ++i) {
    HIC_CHECK(hipStreamWaitEvent(pool.streams[i], pool.fork_event, 0));
  }

  // Loop over zonal wavenumbers
  for (int m = 0; m < n_modes; ++m) {
    hipStream_t mstream = pool.streams[m % pool.n];
    hipblasHandle_t handle = pool.handles[m % pool.n];

    int nlevels = levels[m];
    int lbeta = betalen_max[m];
    const Real *Am = A + a_offsets[m];
    Real *Cm = C + c_offsets[m];
    Real *mwork = work + (m % pool.n) * work_len;
    Real *beta_buf[2] = {mwork, mwork + (size_t)lbeta * n_flds};
    Real *vec_base = mwork + (size_t)2 * lbeta * n_flds;

    // Flat index of NODE(j,k) on level lev of this mode (j and k are 1-based, as in Fortran)
    auto node_index = [&](int lev, int j, int k) {
      int l = lev_offset[m] + lev;
      return lev_node_offset[l] + (k - 1) * lev_ij[l] + (j - 1);
    };

    if (is_transposed) {
      for (int lev = nlevels; lev >= 0; --lev) {
        int l = lev_offset[m] + lev;
        Real *beta = beta_buf[lev % 2];
        Real *beta_prev = beta_buf[(lev + 1) % 2]; // Buffer for level lev-1
        size_t vec_off = 0;
        for (int j = 1; j <= lev_ij[l]; ++j) {
          for (int k = 1; k <= lev_ik[l]; ++k) {
            int idx = node_index(lev, j, k);
            int btst = lev_node_ioffbeta[idx];
            int rank = lev_node_irank[idx];
            int icols = lev_node_icols[idx];
            int n = icols - rank;
            const int *iclist = lev_node_iclist + lev_node_iclist_offset[idx];
            Real *vec = vec_base + vec_off;
            int ld_vec = icols;
            vec_off += (size_t)icols * n_flds;
            check_rank(rank);

            if (lev > 0 && lev == nlevels) {
              // ZBETA = B^T * PVECIN(IFR:ILR,:), where PVECIN is stored transposed
              int irows = lev_node_irows[idx];
              gemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_T, rank, n_flds, irows, (Real)1.0,
                   lev_node_b + lev_node_b_offset[idx], irows,
                   Am + (size_t)(lev_node_ifrow[idx] - 1) * lda, lda, (Real)0.0,
                   beta + btst, lbeta);
            }

            if (n > 0) {
              // ZVECOUT(IRANK+1:ICOLS,:) = PNONIM^T * ZBETA
              gemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_N, n, n_flds, rank, (Real)1.0,
                   lev_node_pnonim + lev_node_pnonim_offset[idx], rank, beta + btst, lbeta,
                   (Real)0.0, vec + rank, ld_vec);
            }

            if (lev == 0) {
              // Scatter into PVECOUT(IFR:ILR,:), stored transposed
              scatter_kernel<<<num_blocks(icols * n_flds), block_size, 0, mstream>>>(
                icols, rank, n_flds, iclist, beta + btst, lbeta, vec, ld_vec,
                Cm + (size_t)(lev_node_ifcol[idx] - 1) * ldc, ldc, 1,
                0, icols, 0, icols, false);
            } else {
              // Scatter into the beta slots of the two children on level lev-1. Odd j assigns,
              // even j (the partner sharing the same children) accumulates.
              int jc = (j + 1) / 2;
              int il = node_index(lev - 1, jc, 2 * k - 1);
              int irankl = lev_node_irank[il];
              int btstl = lev_node_ioffbeta[il];
              int irankr = 0, btstr = 0;
              if (2 * k <= lev_ik[l - 1]) {
                int ir = node_index(lev - 1, jc, 2 * k);
                irankr = lev_node_irank[ir];
                btstr = lev_node_ioffbeta[ir];
              }
              scatter_kernel<<<num_blocks(icols * n_flds), block_size, 0, mstream>>>(
                icols, rank, n_flds, iclist, beta + btst, lbeta, vec, ld_vec,
                beta_prev, 1, lbeta,
                btstl, irankl, btstr, irankl + irankr, j % 2 == 0);
            }
            HIC_CHECK(hipGetLastError());
          }
        }
      }
    } else {
      for (int lev = 0; lev <= nlevels; ++lev) {
        int l = lev_offset[m] + lev;
        Real *beta = beta_buf[lev % 2];
        Real *beta_prev = beta_buf[(lev + 1) % 2]; // Buffer for level lev-1
        size_t vec_off = 0;
        for (int j = 1; j <= lev_ij[l]; ++j) {
          for (int k = 1; k <= lev_ik[l]; ++k) {
            int idx = node_index(lev, j, k);
            int btst = lev_node_ioffbeta[idx];
            int rank = lev_node_irank[idx];
            int icols = lev_node_icols[idx];
            int n = icols - rank;
            const int *iclist = lev_node_iclist + lev_node_iclist_offset[idx];
            Real *vec = vec_base + vec_off;
            int ld_vec = icols;
            vec_off += (size_t)icols * n_flds;
            check_rank(rank);

            if (lev == 0) {
              // Gather from PVECIN(IFR:,:), stored transposed
              gather_kernel<<<num_blocks(icols * n_flds), block_size, 0, mstream>>>(
                icols, rank, n_flds, iclist,
                Am + (size_t)(lev_node_ifcol[idx] - 1) * lda, lda, 1,
                beta + btst, lbeta, vec, ld_vec);
            } else {
              // Gather from the beta slots of the children on level lev-1, which are contiguous
              // starting from the left child
              int il = node_index(lev - 1, (j + 1) / 2, 2 * k - 1);
              gather_kernel<<<num_blocks(icols * n_flds), block_size, 0, mstream>>>(
                icols, rank, n_flds, iclist,
                beta_prev + lev_node_ioffbeta[il], 1, lbeta,
                beta + btst, lbeta, vec, ld_vec);
            }
            HIC_CHECK(hipGetLastError());

            if (n > 0) {
              // ZBETA += PNONIM * ZVECIN(IRANK+1:ICOLS,:)
              gemm(handle, HIPBLAS_OP_N, HIPBLAS_OP_N, rank, n_flds, n, (Real)1.0,
                   lev_node_pnonim + lev_node_pnonim_offset[idx], rank, vec + rank, ld_vec,
                   (Real)1.0, beta + btst, lbeta);
            }

            if (lev == nlevels) {
              // PVECOUT(IFR:ILR,:) = B * ZBETA, where PVECOUT is stored transposed, i.e.
              // PVECOUT^T = ZBETA^T * B^T
              int irows = lev_node_irows[idx];
              gemm(handle, HIPBLAS_OP_T, HIPBLAS_OP_T, n_flds, irows, rank, (Real)1.0,
                   beta + btst, lbeta, lev_node_b + lev_node_b_offset[idx], irows, (Real)0.0,
                   Cm + (size_t)(lev_node_ifrow[idx] - 1) * ldc, ldc);
            }
          }
        }
      }
    }
  }

  // Join the worker streams back onto the caller's stream
  for (int i = 0; i < pool.n; ++i) {
    HIC_CHECK(hipEventRecord(pool.join_events[i], pool.streams[i]));
    HIC_CHECK(hipStreamWaitEvent(stream, pool.join_events[i], 0));
  }
}

// -------------------------------------------------------------------------------------------------
// Fortran-callable wrappers for the mult_butm function
// -------------------------------------------------------------------------------------------------

extern "C" {
void mult_butm_sp(
  char transpose, int n_modes, int n_flds, const int *order, const int *levels,
  const int *betalen_max, const int *lev_offset, const int *lev_ij, const int *lev_ik,
  const int *lev_ibetalen, const int *lev_node_offset, const int *lev_node_ifcol,
  const int *lev_node_ilcol, const int *lev_node_ifrow, const int *lev_node_ilrow,
  const int *lev_node_icols, const int *lev_node_irows, const int *lev_node_irank,
  const int *lev_node_ioffbeta, const int *lev_node_iclist_offset, const int *lev_node_iclist,
  const int *lev_node_pnonim_offset, const float *lev_node_pnonim, const int *lev_node_b_offset,
  const float *lev_node_b, const float *A, int lda, const int64_t *a_offsets, float *C, int ldc,
  const int64_t *c_offsets, const long *stream) {

  mult_butm<float>(transpose, n_modes, n_flds, order, levels, betalen_max,
    lev_offset, lev_ij, lev_ik, lev_ibetalen,
    lev_node_offset, lev_node_ifcol, lev_node_ilcol,
    lev_node_ifrow, lev_node_ilrow, lev_node_icols,
    lev_node_irows, lev_node_irank, lev_node_ioffbeta,
    lev_node_iclist_offset, lev_node_iclist, lev_node_pnonim_offset,
    lev_node_pnonim, lev_node_b_offset, lev_node_b,
    A, lda, a_offsets, C, ldc, c_offsets, reinterpret_cast<hipStream_t>(*stream));
}

void mult_butm_dp(
  char transpose, int n_modes, int n_flds, const int *order, const int *levels,
  const int *betalen_max, const int *lev_offset, const int *lev_ij, const int *lev_ik,
  const int *lev_ibetalen, const int *lev_node_offset, const int *lev_node_ifcol,
  const int *lev_node_ilcol, const int *lev_node_ifrow, const int *lev_node_ilrow,
  const int *lev_node_icols, const int *lev_node_irows, const int *lev_node_irank,
  const int *lev_node_ioffbeta, const int *lev_node_iclist_offset, const int *lev_node_iclist,
  const int *lev_node_pnonim_offset, const double *lev_node_pnonim, const int *lev_node_b_offset,
  const double *lev_node_b, const double *A, int lda, const int64_t *a_offsets, double *C, int ldc,
  const int64_t *c_offsets, const long *stream) {

  mult_butm<double>(transpose, n_modes, n_flds, order, levels, betalen_max,
    lev_offset, lev_ij, lev_ik, lev_ibetalen,
    lev_node_offset, lev_node_ifcol, lev_node_ilcol,
    lev_node_ifrow, lev_node_ilrow, lev_node_icols,
    lev_node_irows, lev_node_irank, lev_node_ioffbeta,
    lev_node_iclist_offset, lev_node_iclist, lev_node_pnonim_offset,
    lev_node_pnonim, lev_node_b_offset, lev_node_b,
    A, lda, a_offsets, C, ldc, c_offsets, reinterpret_cast<hipStream_t>(*stream));
}
}
