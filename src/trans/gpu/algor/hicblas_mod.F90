! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE HICBLAS_MOD

IMPLICIT NONE

INTERFACE
  SUBROUTINE HIP_DGEMM_BATCHED( &
    & CTA, CTB,                 &
    & M, N, K,                  &
    & ALPHA,                    &
    & A, LDA, TDA,              &
    & B, LDB, TDB,              &
    & BETA,                     &
    & C, LDC, TDC,              &
    & BATCHCOUNT, STREAM, ALLOC &
  &) BIND(C, NAME='hipblas_dgemm_wrapper')
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE, C_SIZE_T, C_PTR
    CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
    INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
    REAL(C_DOUBLE),      VALUE            :: ALPHA,BETA
    REAL(C_DOUBLE),      DIMENSION(LDA,*) :: A
    REAL(C_DOUBLE),      DIMENSION(LDB,*) :: B
    REAL(C_DOUBLE),      DIMENSION(LDC,*) :: C
    INTEGER(KIND=C_SIZE_T) :: STREAM
    TYPE(C_PTR), INTENT(IN), VALUE :: ALLOC
  END SUBROUTINE HIP_DGEMM_BATCHED

  SUBROUTINE HIP_SGEMM_BATCHED( &
    & CTA, CTB,                 &
    & M, N, K,                  &
    & ALPHA,                    &
    & A, LDA, TDA,              &
    & B, LDB, TDB,              &
    & BETA,                     &
    & C, LDC, TDC,              &
    & BATCHCOUNT, STREAM, ALLOC &
  &) BIND(C, NAME='hipblas_sgemm_wrapper')
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_FLOAT, C_SIZE_T, C_PTR
    CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
    INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
    REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
    REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
    REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
    REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    INTEGER(KIND=C_SIZE_T) :: STREAM
    TYPE(C_PTR), INTENT(IN), VALUE :: ALLOC
  END SUBROUTINE HIP_SGEMM_BATCHED
END INTERFACE
INTERFACE
  SUBROUTINE CLEAN_GEMM(RESOL_ID) BIND(C, NAME="clean_gemm")
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(IN), VALUE :: RESOL_ID
  END SUBROUTINE
END INTERFACE

INTERFACE
  SUBROUTINE HIP_DGEMM_GROUPED( &
    & RESOL_ID, BLAS_ID, CTA, CTB,        &
    & M, N, K,                  &
    & ALPHA,                    &
    & A, LDA, OFFSETA,          &
    & B, LDB, OFFSETB,          &
    & BETA,                     &
    & C, LDC, OFFSETC,          &
    & BATCHCOUNT, STREAM, ALLOC &
  &) BIND(C, NAME='hipblas_dgemm_wrapper_grouped')
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE, C_SIZE_T, C_PTR, C_INT64_T
    CHARACTER(1,C_CHAR), VALUE :: CTA, CTB
    INTEGER(C_INT), VALUE  :: RESOL_ID, BLAS_ID, M, LDA, LDC, BATCHCOUNT
    INTEGER(C_INT)         :: N(*), K(*), LDB(*)
    INTEGER(C_INT64_T)     :: OFFSETA(*), OFFSETB(*), OFFSETC(*)
    REAL(C_DOUBLE), VALUE  :: ALPHA,BETA
    REAL(C_DOUBLE)         :: A(*), B(*), C(*)
    INTEGER(KIND=C_SIZE_T) :: STREAM
    TYPE(C_PTR), INTENT(IN), VALUE :: ALLOC
  END SUBROUTINE HIP_DGEMM_GROUPED

  SUBROUTINE HIP_SGEMM_GROUPED( &
    & RESOL_ID, BLAS_ID, CTA, CTB,        &
    & M, N, K,                  &
    & ALPHA,                    &
    & A, LDA, OFFSETA,          &
    & B, LDB, OFFSETB,          &
    & BETA,                     &
    & C, LDC, OFFSETC,          &
    & BATCHCOUNT, STREAM, ALLOC &
  &) BIND(C, NAME='hipblas_sgemm_wrapper_grouped')
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_FLOAT, C_SIZE_T, C_PTR, C_INT64_T
    CHARACTER(1,C_CHAR), VALUE :: CTA, CTB
    INTEGER(C_INT), VALUE :: RESOL_ID, BLAS_ID, M, LDA, LDC, BATCHCOUNT
    INTEGER(C_INT)        :: N(*), K(*), LDB(*)
    INTEGER(C_INT64_T)    :: OFFSETA(*), OFFSETB(*), OFFSETC(*)
    REAL(C_FLOAT), VALUE  :: ALPHA,BETA
    REAL(C_FLOAT)         :: A(*), B(*), C(*)
    INTEGER(KIND=C_SIZE_T) :: STREAM
    TYPE(C_PTR), INTENT(IN), VALUE :: ALLOC
  END SUBROUTINE HIP_SGEMM_GROUPED
END INTERFACE

END MODULE HICBLAS_MOD
