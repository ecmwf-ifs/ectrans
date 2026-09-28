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
// layout (row fastest, leading dimensions IBETALEN_MAX and N_ORDER respectively) so the logic
// below maps directly onto the Fortran.
//
// All index/offset arrays are host arrays, apart from lev_node_iclist, lev_node_pnonim and
// lev_node_b which, together with A and C, are device arrays. Indices stored in the flattened
// struct (IFCOL, IFROW, ICLIST) are Fortran 1-based; IOFFBETA and all *_OFFSET arrays are 0-based.

#include <algorithm>
#include <cstdint>
#include <cstdlib>

#include "hicblas.h"

namespace {

hipblasHandle_t get_mult_butm_handle() {
  static hipblasHandle_t handle;
  static bool initialised = false;
  if (!initialised) {
    HICBLAS_CHECK(hipblasCreate(&handle));
    initialised = true;
  }
  return handle;
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

  hipblasHandle_t handle = get_mult_butm_handle();
  HICBLAS_CHECK(hipblasSetStream(handle, stream));

  // Work space: two "beta" buffers (ZBETA(:,:,0:1)) and one vector buffer (ZVECIN in the normal
  // case, ZVECOUT in the transposed case), sized for the largest mode
  int lbeta_max = 0, order_max = 0;
  for (int m = 0; m < n_modes; ++m) {
    lbeta_max = std::max(lbeta_max, betalen_max[m]);
    order_max = std::max(order_max, order[m]);
  }
  Real *work = get_workspace<Real>((size_t)(2 * lbeta_max + order_max) * n_flds);
  Real *beta_buf[2] = {work, work + (size_t)lbeta_max * n_flds};
  Real *vec = work + (size_t)2 * lbeta_max * n_flds;

  // Loop over zonal wavenumbers
  for (int m = 0; m < n_modes; ++m) {
    int nlevels = levels[m];
    int lbeta = betalen_max[m];
    int ld_vec = order[m];
    const Real *Am = A + a_offsets[m];
    Real *Cm = C + c_offsets[m];

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
        for (int j = 1; j <= lev_ij[l]; ++j) {
          for (int k = 1; k <= lev_ik[l]; ++k) {
            int idx = node_index(lev, j, k);
            int btst = lev_node_ioffbeta[idx];
            int rank = lev_node_irank[idx];
            int icols = lev_node_icols[idx];
            int n = icols - rank;
            const int *iclist = lev_node_iclist + lev_node_iclist_offset[idx];
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
              scatter_kernel<<<num_blocks(icols * n_flds), block_size, 0, stream>>>(
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
              scatter_kernel<<<num_blocks(icols * n_flds), block_size, 0, stream>>>(
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
        for (int j = 1; j <= lev_ij[l]; ++j) {
          for (int k = 1; k <= lev_ik[l]; ++k) {
            int idx = node_index(lev, j, k);
            int btst = lev_node_ioffbeta[idx];
            int rank = lev_node_irank[idx];
            int icols = lev_node_icols[idx];
            int n = icols - rank;
            const int *iclist = lev_node_iclist + lev_node_iclist_offset[idx];
            check_rank(rank);

            if (lev == 0) {
              // Gather from PVECIN(IFR:,:), stored transposed
              gather_kernel<<<num_blocks(icols * n_flds), block_size, 0, stream>>>(
                icols, rank, n_flds, iclist,
                Am + (size_t)(lev_node_ifcol[idx] - 1) * lda, lda, 1,
                beta + btst, lbeta, vec, ld_vec);
            } else {
              // Gather from the beta slots of the children on level lev-1, which are contiguous
              // starting from the left child
              int il = node_index(lev - 1, (j + 1) / 2, 2 * k - 1);
              gather_kernel<<<num_blocks(icols * n_flds), block_size, 0, stream>>>(
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
