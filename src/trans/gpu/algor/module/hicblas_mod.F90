! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

#ifdef HIPGPU
#define hipblasSgemm 'hipblasSgemm'
#define hipblasDgemm 'hipblasDgemm'
#elif defined CUDAGPU
#define hipblasSgemm 'cublasSgemm'
#define hipblasDgemm 'cublasDgemm'
#endif


MODULE HICBLAS_MOD

USE PARKIND1, ONLY : JPIM, JPRM, JPRD
USE ISO_C_BINDING

!
! Define the interfaces to HIP/CUDA C code via a common wrapper interface
!
interface hip_gemm
!
! void hipblasSgemm (char transa, char transb, int m, int n,
! int k, float alpha, const float *A, int lda,
! const float *B, int ldb, float beta, float *C, int ldc)
!
subroutine hip_sgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='hipblasSgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_float),value :: alpha,beta
real(c_float), dimension(lda,*) :: A
real(c_float), dimension(ldb,*) :: B
real(c_float), dimension(ldc,*) :: C
end subroutine hip_sgemm

!
! void hipblasDgemm (char transa, char transb, int m, int n,
! int k, double alpha, const double *A, int lda,
! const double *B, int ldb, double beta, double *C, int ldc)
!
subroutine hip_dgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='hipblasDgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_double),value :: alpha,beta
real(c_double), dimension(lda,*) :: A
real(c_double), dimension(ldb,*) :: B
real(c_double), dimension(ldc,*) :: C
end subroutine hip_dgemm
end interface


INTERFACE
    SUBROUTINE HIP_DGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasDgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_DOUBLE),      VALUE            :: ALPHA,BETA
        REAL(C_DOUBLE),      DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),      DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),      DIMENSION(LDC,*) :: C
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
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasDgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_INT),       VALUE            :: TDA,TDB,TDC
        REAL(C_DOUBLE),       VALUE            :: ALPHA, BETA
        REAL(C_DOUBLE),       DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),       DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),       DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_DGEMM_STRIDED_BATCHED
END INTERFACE

INTERFACE

    subroutine hip_dgemm_batched_finalize() bind(C,name='hipblasDgemmBatched_finalize')
    end subroutine hip_dgemm_batched_finalize

END INTERFACE

INTERFACE

    SUBROUTINE HIP_SGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasSgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
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
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasSgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_INT),       VALUE            :: TDA,TDB,TDC
        REAL(C_FLOAT),        VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),        DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),        DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),        DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_SGEMM_STRIDED_BATCHED

END INTERFACE

INTERFACE
    subroutine hip_sgemm_batched_finalize() bind(C,name='hipblasSgemmBatched_finalize')
    end subroutine hip_sgemm_batched_finalize


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
END INTERFACE

INTERFACE
SUBROUTINE HIP_DGEMM_GROUPED(&
    & CTA, CTB,               &
    & M, N, K,                &
    & ALPHA,                  &
    & A, LDA, TDA,            &
    & B, LDB, TDB,            &
    & BETA,                   &
    & C, LDC, TDC,            &
    & BATCHCOUNT              &
&) BIND(C, NAME='blas_dgemm_wrapper_grouped')
    USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
!!    INTEGER(C_INT), VALUE :: CTA, CTB, M, N(:), K(:), LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
    INTEGER(C_INT), VALUE :: M, N(:), K(:), LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
    REAL(C_DOUBLE), VALUE  :: ALPHA,BETA
    REAL(C_DOUBLE)         :: A(*), B(*), C(*)
END SUBROUTINE HIP_DGEMM_GROUPED
SUBROUTINE HIP_SGEMM_GROUPED(&
    & CTA, CTB,               &
    & M, N, K,                &
    & ALPHA,                  &
    & A, LDA, TDA,            &
    & B, LDB, TDB,            &
    & BETA,                   &
    & C, LDC, TDC,            &
    & BATCHCOUNT              &
&) BIND(C, NAME='blas_sgemm_wrapper_grouped')
    USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
!!    INTEGER(C_INT), VALUE :: CTA, CTB, M, N(:), K(:), LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
    INTEGER(C_INT), VALUE :: M, N(:), K(:), LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
    REAL(C_FLOAT), VALUE  :: ALPHA,BETA
    REAL(C_FLOAT)         :: A(*), B(*), C(*)
END SUBROUTINE HIP_SGEMM_GROUPED
END INTERFACE

CONTAINS

SUBROUTINE HIP_DGEMM_GROUPED_OVERLOAD(&
    & TRANSA, TRANSB, &
    & M, N, K, &
    & ALPHA, &
    & AARRAY, LDA, STRIDEA, &
    & BARRAY, LDB, STRIDEB, &
    & BETA, &
    & CARRAY, LDC, STRIDEC, &
    & BATCHCOUNT)
!!  INTEGER(KIND=C_INT), INTENT(IN) :: TRANSA
!!  INTEGER(KIND=C_INT), INTENT(IN) :: TRANSB
  CHARACTER(1,C_CHAR), VALUE :: TRANSA, TRANSB
  INTEGER(KIND=JPIM) :: M
  INTEGER(KIND=JPIM) :: N(:)
  INTEGER(KIND=JPIM) :: K(:)
  REAL(KIND=JPRD) :: ALPHA
  REAL(KIND=JPRD), DIMENSION(*) :: AARRAY
  INTEGER(KIND=JPIM) :: LDA
  INTEGER(KIND=JPIM) :: STRIDEA
  REAL(KIND=JPRD), DIMENSION(*) :: BARRAY
  INTEGER(KIND=JPIM) :: LDB
  INTEGER(KIND=JPIM) :: STRIDEB
  REAL(KIND=JPRD) :: BETA
  REAL(KIND=JPRD), DIMENSION(*) :: CARRAY
  INTEGER(KIND=JPIM) :: LDC
  INTEGER(KIND=JPIM) :: STRIDEC
  INTEGER(KIND=JPIM) :: BATCHCOUNT

  CALL HIP_DGEMM_GROUPED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)

END SUBROUTINE HIP_DGEMM_GROUPED_OVERLOAD

SUBROUTINE HIP_SGEMM_GROUPED_OVERLOAD(&
    & TRANSA, TRANSB, & 
    & M, N, K, & 
    & ALPHA, & 
    & AARRAY, LDA, STRIDEA, & 
    & BARRAY, LDB, STRIDEB, & 
    & BETA, & 
    & CARRAY, LDC, STRIDEC, & 
    & BATCHCOUNT)
!!  INTEGER(KIND=C_INT), INTENT(IN) :: TRANSA
!!  INTEGER(KIND=C_INT), INTENT(IN) :: TRANSB
  CHARACTER(1,C_CHAR), VALUE :: TRANSA, TRANSB
  INTEGER(KIND=JPIM) :: M 
  INTEGER(KIND=JPIM) :: N(:)
  INTEGER(KIND=JPIM) :: K(:)
  REAL(KIND=JPRM) :: ALPHA
  REAL(KIND=JPRM), DIMENSION(*) :: AARRAY
  INTEGER(KIND=JPIM) :: LDA
  INTEGER(KIND=JPIM) :: STRIDEA
  REAL(KIND=JPRM), DIMENSION(*) :: BARRAY
  INTEGER(KIND=JPIM) :: LDB
  INTEGER(KIND=JPIM) :: STRIDEB
  REAL(KIND=JPRM) :: BETA
  REAL(KIND=JPRM), DIMENSION(*) :: CARRAY
  INTEGER(KIND=JPIM) :: LDC
  INTEGER(KIND=JPIM) :: STRIDEC
  INTEGER(KIND=JPIM) :: BATCHCOUNT
    
  CALL HIP_SGEMM_GROUPED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)

END SUBROUTINE HIP_SGEMM_GROUPED_OVERLOAD

END MODULE HICBLAS_MOD
