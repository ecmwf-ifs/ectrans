! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

#if defined CUDAGPU
#define hipblasSgemm 'cublasSgemm'
#define hipblasDgemm 'cublasDgemm'
#define ACC_GET_HIP_STREAM ACC_GET_CUDA_STREAM
#define OPENACC_LIB OPENACC
#endif

MODULE HICBLAS_MOD

USE PARKIND1, ONLY : JPIM, JPRM, JPRD
USE GROWING_ALLOCATOR_MOD, ONLY: GROWING_ALLOCATION_TYPE
USE ISO_C_BINDING
USE OPENACC_LIB, ONLY: ACC_GET_HIP_STREAM

IMPLICIT NONE

  INTERFACE HIP_GEMM_BATCHED
    MODULE PROCEDURE HIP_DGEMM_BATCHED_OVERLOAD
    MODULE PROCEDURE HIP_SGEMM_BATCHED_OVERLOAD
    MODULE PROCEDURE HIP_DGEMM_GROUPED_OVERLOAD
    MODULE PROCEDURE HIP_SGEMM_GROUPED_OVERLOAD
  END INTERFACE HIP_GEMM_BATCHED

!
! Define the interfaces to HIP/CUDA C code via a common wrapper interface
!
interface hip_gemm
!
! void hipblasSgemm (char transa, char transb, int m, int n,
! int k, float alpha, const float *A, int lda,
! const float *B, int ldb, float beta, float *C, int ldc)
!
SUBROUTINE HIP_SGEMM(CTA, CTB, M, N, K,&
ALPHA, A, LDA, B, LDB, BETA, C, LDC) BIND(C,NAME='hipblasSgemm')
USE ISO_C_BINDING
CHARACTER(1,C_CHAR),VALUE :: CTA, CTB
INTEGER(C_INT),     VALUE :: M,N,K,LDA,LDB,LDC
REAL(C_FLOAT),      VALUE :: ALPHA,BETA
REAL(C_FLOAT), DIMENSION(LDA,*) :: A
REAL(C_FLOAT), DIMENSION(LDB,*) :: B
REAL(C_FLOAT), DIMENSION(LDC,*) :: C
END SUBROUTINE HIP_SGEMM

!
! void hipblasDgemm (char transa, char transb, int m, int n,
! int k, double alpha, const double *A, int lda,
! const double *B, int ldb, double beta, double *C, int ldc)
!
SUBROUTINE HIP_DGEMM(CTA, CTB, M, N, K,&
ALPHA, A, LDA, B, LDB, BETA, C, LDC) BIND(C,NAME='hipblasDgemm')
USE ISO_C_BINDING
CHARACTER(1,C_CHAR),VALUE :: CTA, CTB
INTEGER(C_INT),     VALUE :: M,N,K,LDA,LDB,LDC
REAL(C_DOUBLE),     VALUE :: ALPHA,BETA
REAL(C_DOUBLE), DIMENSION(LDA,*) :: A
REAL(C_DOUBLE), DIMENSION(LDB,*) :: B
REAL(C_DOUBLE), DIMENSION(LDC,*) :: C
END SUBROUTINE HIP_DGEMM
END INTERFACE

INTERFACE
    SUBROUTINE HIP_DGEMM_BATCHED(   &
        & CTA, CTB,                 &
        & M, N, K,                  &
        & ALPHA,                    &
        & A, LDA, TDA,              &
        & B, LDB, TDB,              &
        & BETA,                     &
        & C, LDC, TDC,              &
        & BATCHCOUNT, STREAM, ALLOC &
    &) BIND(C, NAME='hipblas_dgemm_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_DOUBLE),      VALUE            :: ALPHA,BETA
        REAL(C_DOUBLE),      DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),      DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),      DIMENSION(LDC,*) :: C
        INTEGER(KIND=C_SIZE_T) :: STREAM
        TYPE(C_PTR), INTENT(IN), VALUE :: ALLOC
    END SUBROUTINE HIP_DGEMM_BATCHED
END INTERFACE

INTERFACE
    SUBROUTINE HIP_DGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT, STREAM      &
    &) BIND(C, NAME='hipblasDgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_INT),       VALUE            :: TDA,TDB,TDC
        REAL(C_DOUBLE),       VALUE            :: ALPHA, BETA
        REAL(C_DOUBLE),       DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),       DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),       DIMENSION(LDC,*) :: C
        INTEGER(KIND=C_SIZE_T) :: STREAM
    END SUBROUTINE HIP_DGEMM_STRIDED_BATCHED
END INTERFACE

INTERFACE
    SUBROUTINE HIP_DGEMM_BATCHED_FINALIZE() BIND(C,NAME='hipblasDgemmBatched_finalize')
    END SUBROUTINE HIP_DGEMM_BATCHED_FINALIZE
END INTERFACE

INTERFACE
    SUBROUTINE HIP_SGEMM_BATCHED(   &
        & CTA, CTB,                 &
        & M, N, K,                  &
        & ALPHA,                    &
        & A, LDA, TDA,              &
        & B, LDB, TDB,              &
        & BETA,                     &
        & C, LDC, TDC,              &
        & BATCHCOUNT, STREAM, ALLOC &
    &) BIND(C, NAME='hipblas_sgemm_wrapper')
        USE ISO_C_BINDING
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
    SUBROUTINE HIP_SGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT, STREAM      &
    &) BIND(C, NAME='hipblasSgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_INT),       VALUE            :: TDA,TDB,TDC
        REAL(C_FLOAT),        VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),        DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),        DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),        DIMENSION(LDC,*) :: C
        INTEGER(KIND=C_SIZE_T) :: STREAM
    END SUBROUTINE HIP_SGEMM_STRIDED_BATCHED
END INTERFACE

INTERFACE
    SUBROUTINE HIP_SGEMM_BATCHED_FINALIZE() BIND(C,NAME='hipblasSgemmBatched_finalize')
    END SUBROUTINE HIP_SGEMM_BATCHED_FINALIZE
END INTERFACE

INTERFACE
    SUBROUTINE HIP_STCGEMM_BATCHED(&
        & CTA, CTB,                &
        & M, N, K,                 &
        & ALPHA,                   &
        & A, LDA, TDA,             &
        & B, LDB, TDB,             &
        & BETA,                    &
        & C, LDC, TDC,             &
        & BATCHCOUNT               &
    &) BIND(C, NAME='hipblasSTCgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_STCGEMM_BATCHED
END INTERFACE

INTERFACE
SUBROUTINE HIP_DGEMM_GROUPED(   &
    & BLAS_ID, CTA, CTB,        &
    & M, N, K,                  &
    & ALPHA,                    &
    & A, LDA, OFFSETA,          &
    & B, LDB, OFFSETB,          &
    & BETA,                     &
    & C, LDC, OFFSETC,          &
    & BATCHCOUNT, STREAM, ALLOC &
&) BIND(C, NAME='blas_dgemm_wrapper_grouped')
    USE ISO_C_BINDING
    CHARACTER(1,C_CHAR), VALUE :: CTA, CTB
    INTEGER(C_INT), VALUE  :: BLAS_ID, M, LDA, LDB, LDC, BATCHCOUNT
    INTEGER(C_INT)         :: N(*), K(*), OFFSETA(*), OFFSETB(*), OFFSETC(*)
    REAL(C_DOUBLE), VALUE  :: ALPHA,BETA
    REAL(C_DOUBLE)         :: A(*), B(*), C(*)
    INTEGER(KIND=C_SIZE_T) :: STREAM
    TYPE(C_PTR), INTENT(IN), VALUE :: ALLOC
END SUBROUTINE HIP_DGEMM_GROUPED
SUBROUTINE HIP_SGEMM_GROUPED(   &
    & BLAS_ID, CTA, CTB,        &
    & M, N, K,                  &
    & ALPHA,                    &
    & A, LDA, OFFSETA,          &
    & B, LDB, OFFSETB,          &
    & BETA,                     &
    & C, LDC, OFFSETC,          &
    & BATCHCOUNT, STREAM, ALLOC &
&) BIND(C, NAME='blas_sgemm_wrapper_grouped')
    USE ISO_C_BINDING
    CHARACTER(1,C_CHAR), VALUE :: CTA, CTB
    INTEGER(C_INT), VALUE :: BLAS_ID, M, LDA, LDB, LDC, BATCHCOUNT
    INTEGER(C_INT)        :: N(*), K(*), OFFSETA(*), OFFSETB(*), OFFSETC(*)
    REAL(C_FLOAT), VALUE  :: ALPHA,BETA
    REAL(C_FLOAT)         :: A(*), B(*), C(*)
    INTEGER(KIND=C_SIZE_T) :: STREAM
    TYPE(C_PTR), INTENT(IN), VALUE :: ALLOC
END SUBROUTINE HIP_SGEMM_GROUPED
END INTERFACE

CONTAINS

SUBROUTINE HIP_DGEMM_BATCHED_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT, STREAM, ALLOC)
    CHARACTER(1,C_CHAR), VALUE :: TRANSA, TRANSB
    INTEGER(KIND=JPIM) :: M
    INTEGER(KIND=JPIM) :: N
    INTEGER(KIND=JPIM) :: K
    REAL(KIND=JPRD) :: ALPHA
    REAL(KIND=JPRD), DIMENSION(:) :: AARRAY
    INTEGER(KIND=JPIM) :: LDA
    INTEGER(KIND=JPIM) :: STRIDEA
    REAL(KIND=JPRD), DIMENSION(:,:) :: BARRAY
    INTEGER(KIND=JPIM) :: LDB
    INTEGER(KIND=JPIM) :: STRIDEB
    REAL(KIND=JPRD) :: BETA
    REAL(KIND=JPRD), DIMENSION(:) :: CARRAY
    INTEGER(KIND=JPIM) :: LDC
    INTEGER(KIND=JPIM) :: STRIDEC
    INTEGER(KIND=JPIM) :: BATCHCOUNT
    INTEGER(KIND=C_INT) :: STREAM
    TYPE(GROWING_ALLOCATION_TYPE), INTENT(IN) :: ALLOC

    INTEGER(KIND=C_LONG) :: HIP_STREAM

    HIP_STREAM = INT(ACC_GET_HIP_STREAM(STREAM), C_LONG)

#if defined(_CRAYFTN)
    !$ACC HOST_DATA USE_DEVICE(AARRAY,BARRAY,CARRAY)
#endif
    CALL HIP_DGEMM_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT, HIP_STREAM, C_LOC(ALLOC))
#if defined(_CRAYFTN)
    !$ACC END HOST_DATA
#endif
  END SUBROUTINE HIP_DGEMM_BATCHED_OVERLOAD

  SUBROUTINE HIP_SGEMM_BATCHED_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT, STREAM, ALLOC)
    CHARACTER(1,C_CHAR), VALUE :: TRANSA, TRANSB
    INTEGER(KIND=JPIM) :: M
    INTEGER(KIND=JPIM) :: N
    INTEGER(KIND=JPIM) :: K
    REAL(KIND=JPRM) :: ALPHA
    REAL(KIND=JPRM), DIMENSION(:) :: AARRAY
    INTEGER(KIND=JPIM) :: LDA
    INTEGER(KIND=JPIM) :: STRIDEA
    REAL(KIND=JPRM), DIMENSION(*) :: BARRAY
    INTEGER(KIND=JPIM) :: LDB
    INTEGER(KIND=JPIM) :: STRIDEB
    REAL(KIND=JPRM) :: BETA
    REAL(KIND=JPRM), DIMENSION(:) :: CARRAY
    INTEGER(KIND=JPIM) :: LDC
    INTEGER(KIND=JPIM) :: STRIDEC
    INTEGER(KIND=JPIM) :: BATCHCOUNT
    INTEGER(KIND=C_INT) :: STREAM
    TYPE(GROWING_ALLOCATION_TYPE), INTENT(IN) :: ALLOC

    INTEGER(KIND=C_LONG) :: HIP_STREAM

    HIP_STREAM = INT(ACC_GET_HIP_STREAM(STREAM), C_LONG)

    CALL HIP_SGEMM_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT, HIP_STREAM, C_LOC(ALLOC))
  END SUBROUTINE HIP_SGEMM_BATCHED_OVERLOAD


  SUBROUTINE HIP_DGEMM_GROUPED_OVERLOAD( &
      & BLAS_ID, TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, OFFSETA, &
      & BARRAY, LDB, OFFSETB, &
      & BETA, &
      & CARRAY, LDC, OFFSETC, &
      & BATCHCOUNT, STREAM, ALLOC)
    INTEGER(KIND=C_INT), INTENT(IN) :: BLAS_ID
    CHARACTER(1,C_CHAR), VALUE :: TRANSA, TRANSB
    INTEGER(KIND=JPIM) :: M
    INTEGER(KIND=JPIM) :: N(:)
    INTEGER(KIND=JPIM) :: K(:)
    REAL(KIND=JPRD) :: ALPHA
    REAL(KIND=JPRD), DIMENSION(:) :: AARRAY
    INTEGER(KIND=JPIM) :: LDA
    INTEGER(KIND=JPIM) :: OFFSETA(:)
    REAL(KIND=JPRD), DIMENSION(*) :: BARRAY
    INTEGER(KIND=JPIM) :: LDB
    INTEGER(KIND=JPIM) :: OFFSETB(:)
    REAL(KIND=JPRD) :: BETA
    REAL(KIND=JPRD), DIMENSION(:) :: CARRAY
    INTEGER(KIND=JPIM) :: LDC
    INTEGER(KIND=JPIM) :: OFFSETC(:)
    INTEGER(KIND=JPIM) :: BATCHCOUNT
    INTEGER(KIND=C_INT) :: STREAM
    TYPE(GROWING_ALLOCATION_TYPE), INTENT(IN) :: ALLOC

    INTEGER(KIND=C_LONG) :: HIP_STREAM

    HIP_STREAM = INT(ACC_GET_HIP_STREAM(STREAM), C_LONG)

    CALL HIP_DGEMM_GROUPED( &
      & BLAS_ID, TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, OFFSETA, &
      & BARRAY, LDB, OFFSETB, &
      & BETA, &
      & CARRAY, LDC, OFFSETC, &
      & BATCHCOUNT, HIP_STREAM, C_LOC(ALLOC))

  END SUBROUTINE HIP_DGEMM_GROUPED_OVERLOAD

  SUBROUTINE HIP_SGEMM_GROUPED_OVERLOAD(&
      & BLAS_ID, TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, OFFSETA, &
      & BARRAY, LDB, OFFSETB, &
      & BETA, &
      & CARRAY, LDC, OFFSETC, &
      & BATCHCOUNT, STREAM, ALLOC)
    INTEGER(KIND=C_INT), INTENT(IN) :: BLAS_ID
    CHARACTER(1,C_CHAR), VALUE :: TRANSA, TRANSB
    INTEGER(KIND=JPIM) :: M
    INTEGER(KIND=JPIM) :: N(:)
    INTEGER(KIND=JPIM) :: K(:)
    REAL(KIND=JPRM) :: ALPHA
    REAL(KIND=JPRM), DIMENSION(:) :: AARRAY
    INTEGER(KIND=JPIM) :: LDA
    INTEGER(KIND=JPIM) :: OFFSETA(:)
    REAL(KIND=JPRM), DIMENSION(:,:,:) :: BARRAY
    INTEGER(KIND=JPIM) :: LDB
    INTEGER(KIND=JPIM) :: OFFSETB(:)
    REAL(KIND=JPRM) :: BETA
    REAL(KIND=JPRM), DIMENSION(:) :: CARRAY
    INTEGER(KIND=JPIM) :: LDC
    INTEGER(KIND=JPIM) :: OFFSETC(:)
    INTEGER(KIND=JPIM) :: BATCHCOUNT
    INTEGER(KIND=C_INT) :: STREAM
    TYPE(GROWING_ALLOCATION_TYPE), INTENT(IN) :: ALLOC

    INTEGER(KIND=C_LONG) :: HIP_STREAM

    HIP_STREAM = INT(ACC_GET_HIP_STREAM(STREAM), C_LONG)

#if defined(_CRAYFTN)
    !$ACC HOST_DATA USE_DEVICE(AARRAY,BARRAY,CARRAY)
#endif
    CALL HIP_SGEMM_GROUPED( &
      & BLAS_ID, TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, OFFSETA, &
      & BARRAY, LDB, OFFSETB, &
      & BETA, &
      & CARRAY, LDC, OFFSETC, &
      & BATCHCOUNT, HIP_STREAM, C_LOC(ALLOC))
#if defined(_CRAYFTN)
    !$ACC END HOST_DATA
#endif

  END SUBROUTINE HIP_SGEMM_GROUPED_OVERLOAD

END MODULE HICBLAS_MOD
