! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE CUDA_GEMM_BATCHED_MOD
  USE CUBLAS_MOD
  USE PARKIND1, ONLY: JPRD, JPRM, JPRL, JPIM, JPIB

  IMPLICIT NONE 

  PRIVATE
  PUBLIC CUDA_GEMM_BATCHED, CUDA_DGEMM_BATCHED !!,CUDA_TCGEMM_BATCHED !! , CUDA_GEMM_STRIDED_BATCHED

  INTERFACE CUDA_GEMM_BATCHED
    MODULE PROCEDURE CUDA_DGEMM_BATCHED_OVERLOAD
    MODULE PROCEDURE CUDA_SGEMM_BATCHED_OVERLOAD
    MODULE PROCEDURE CUDA_DGEMM_STRIDED_BATCHED_OVERLOAD
    MODULE PROCEDURE CUDA_SGEMM_STRIDED_BATCHED_OVERLOAD
    MODULE PROCEDURE CUDA_DGEMM_BATCHED_1D_OVERLOAD
    MODULE PROCEDURE CUDA_SGEMM_BATCHED_1D_OVERLOAD
    MODULE PROCEDURE CUDA_DGEMM_STRIDED_BATCHED_1D_OVERLOAD
    MODULE PROCEDURE CUDA_SGEMM_STRIDED_BATCHED_1D_OVERLOAD
  END INTERFACE CUDA_GEMM_BATCHED

!  INTERFACE CUDA_GEMM_STRIDED_BATCHED
!    MODULE PROCEDURE CUDA_SGEMM_STRIDED_BATCHED_OVERLOAD
!  END INTERFACE CUDA_GEMM_STRIDED_BATCHED

!  INTERFACE CUDA_TCGEMM_BATCHED
!    MODULE PROCEDURE CUDA_DTCGEMM_BATCHED_OVERLOAD
!    MODULE PROCEDURE CUDA_STCGEMM_BATCHED_OVERLOAD
!  END INTERFACE CUDA_TCGEMM_BATCHED

CONTAINS
  SUBROUTINE CUDA_DGEMM_BATCHED_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER,                         INTENT(IN)  :: TRANSA
    CHARACTER,                         INTENT(IN)  :: TRANSB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
    REAL(KIND=JPRD),                   INTENT(IN)  :: ALPHA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(IN)  :: AARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEB
    REAL(KIND=JPRD),                   INTENT(IN)  :: BETA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

     print *, ' dgemm_batched_overload'
     CALL CUDA_DGEMM_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_DGEMM_BATCHED_OVERLOAD

  SUBROUTINE CUDA_SGEMM_BATCHED_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER,                         INTENT(IN)  :: TRANSA
    CHARACTER,                         INTENT(IN)  :: TRANSB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
    REAL(KIND=JPRM),                   INTENT(IN)  :: ALPHA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(IN)  :: AARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEB
    REAL(KIND=JPRM),                   INTENT(IN)  :: BETA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

    CALL CUDA_SGEMM_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_SGEMM_BATCHED_OVERLOAD

  SUBROUTINE CUDA_SGEMM_STRIDED_BATCHED_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER,                         INTENT(IN)  :: TRANSA
    CHARACTER,                         INTENT(IN)  :: TRANSB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
    REAL(KIND=JPRM),                   INTENT(IN)  :: ALPHA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(IN)  :: AARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
    INTEGER(KIND=JPIB),                INTENT(IN)  :: STRIDEA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
    INTEGER(KIND=JPIB),                INTENT(IN)  :: STRIDEB
    REAL(KIND=JPRM),                   INTENT(IN)  :: BETA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
    INTEGER(KIND=JPIB),                INTENT(IN)  :: STRIDEC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

    CALL CUDA_SGEMM_STRIDED_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_SGEMM_STRIDED_BATCHED_OVERLOAD

  SUBROUTINE CUDA_DGEMM_BATCHED_1D_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER,                         INTENT(IN)  :: TRANSA
    CHARACTER,                         INTENT(IN)  :: TRANSB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
    REAL(KIND=JPRD),                   INTENT(IN)  :: ALPHA
    REAL(KIND=JPRD), DIMENSION(:),     INTENT(IN)  :: AARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEB
    REAL(KIND=JPRD),                   INTENT(IN)  :: BETA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

    print *, ' dgemm_batched_1d_overload'
    CALL CUDA_DGEMM_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_DGEMM_BATCHED_1D_OVERLOAD

  SUBROUTINE CUDA_SGEMM_BATCHED_1D_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER,                         INTENT(IN)  :: TRANSA
    CHARACTER,                         INTENT(IN)  :: TRANSB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
    REAL(KIND=JPRM),                   INTENT(IN)  :: ALPHA
    REAL(KIND=JPRM), DIMENSION(:),     INTENT(IN)  :: AARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEB
    REAL(KIND=JPRM),                   INTENT(IN)  :: BETA
    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

    CALL CUDA_SGEMM_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_SGEMM_BATCHED_1D_OVERLOAD

  SUBROUTINE CUDA_SGEMM_STRIDED_BATCHED_1D_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER, INTENT(IN) :: TRANSA
    CHARACTER, INTENT(IN) :: TRANSB
    INTEGER(KIND=JPIM) :: M
    INTEGER(KIND=JPIM) :: N
    INTEGER(KIND=JPIM) :: K
    REAL(KIND=JPRM) :: ALPHA
    REAL(KIND=JPRM), DIMENSION(:) :: AARRAY
    INTEGER(KIND=JPIM) :: LDA
    INTEGER(KIND=JPIB) :: STRIDEA
    REAL(KIND=JPRM), DIMENSION(:,:,:) :: BARRAY
    INTEGER(KIND=JPIM) :: LDB
    INTEGER(KIND=JPIB) :: STRIDEB
    REAL(KIND=JPRM) :: BETA
    REAL(KIND=JPRM), DIMENSION(:) :: CARRAY
    INTEGER(KIND=JPIM) :: LDC
    INTEGER(KIND=JPIB) :: STRIDEC
    INTEGER(KIND=JPIM) :: BATCHCOUNT

    CALL CUDA_SGEMM_STRIDED_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_SGEMM_STRIDED_BATCHED_1D_OVERLOAD


  SUBROUTINE CUDA_DGEMM_STRIDED_BATCHED_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER,                         INTENT(IN)  :: TRANSA
    CHARACTER,                         INTENT(IN)  :: TRANSB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
    REAL(KIND=JPRD),                   INTENT(IN)  :: ALPHA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(IN)  :: AARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
    INTEGER(KIND=JPIB),                INTENT(IN)  :: STRIDEA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
    INTEGER(KIND=JPIB),                INTENT(IN)  :: STRIDEB
    REAL(KIND=JPRD),                   INTENT(IN)  :: BETA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
    INTEGER(KIND=JPIB),                INTENT(IN)  :: STRIDEC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

    print *, ' dgemm_strided_batched_overload'
    CALL CUDA_DGEMM_STRIDED_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_DGEMM_STRIDED_BATCHED_OVERLOAD

  SUBROUTINE CUDA_DGEMM_STRIDED_BATCHED_1D_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER, INTENT(IN) :: TRANSA
    CHARACTER, INTENT(IN) :: TRANSB
    INTEGER(KIND=JPIM) :: M
    INTEGER(KIND=JPIM) :: N
    INTEGER(KIND=JPIM) :: K
    REAL(KIND=JPRD) :: ALPHA
    REAL(KIND=JPRD), DIMENSION(:) :: AARRAY
    INTEGER(KIND=JPIM) :: LDA
    INTEGER(KIND=JPIB) :: STRIDEA
    REAL(KIND=JPRD), DIMENSION(:,:,:) :: BARRAY
    INTEGER(KIND=JPIM) :: LDB
    INTEGER(KIND=JPIB) :: STRIDEB
    REAL(KIND=JPRD) :: BETA
    REAL(KIND=JPRD), DIMENSION(:) :: CARRAY
    INTEGER(KIND=JPIM) :: LDC
    INTEGER(KIND=JPIB) :: STRIDEC
    INTEGER(KIND=JPIM) :: BATCHCOUNT

    print *, ' dgemm_strided_batched_1d_overload'
    CALL CUDA_DGEMM_STRIDED_BATCHED( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
  END SUBROUTINE CUDA_DGEMM_STRIDED_BATCHED_1D_OVERLOAD

  SUBROUTINE CUDA_DTCGEMM_BATCHED_OVERLOAD( &
      & TRANSA, TRANSB, &
      & M, N, K, &
      & ALPHA, &
      & AARRAY, LDA, STRIDEA, &
      & BARRAY, LDB, STRIDEB, &
      & BETA, &
      & CARRAY, LDC, STRIDEC, &
      & BATCHCOUNT)
    CHARACTER,                         INTENT(IN)  :: TRANSA
    CHARACTER,                         INTENT(IN)  :: TRANSB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
    REAL(KIND=JPRD),                   INTENT(IN)  :: ALPHA
    REAL(KIND=JPRL), DIMENSION(:,:,:), INTENT(IN)  :: AARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEA
    REAL(KIND=JPRL), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEB
    REAL(KIND=JPRD),                   INTENT(IN)  :: BETA
    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEC
    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

    STOP "CUDA_TCGEMM_BATCHED not implemented for double-precision"
  END SUBROUTINE CUDA_DTCGEMM_BATCHED_OVERLOAD

!  SUBROUTINE CUDA_STCGEMM_BATCHED_OVERLOAD( &
!      & TRANSA, TRANSB, &
!      & M, N, K, &
!      & ALPHA, &
!      & AARRAY, LDA, STRIDEA, &
!      & BARRAY, LDB, STRIDEB, &
!      & BETA, &
!      & CARRAY, LDC, STRIDEC, &
!      & BATCHCOUNT)
!    CHARACTER,                         INTENT(IN)  :: TRANSA
!    CHARACTER,                         INTENT(IN)  :: TRANSB
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
!    REAL(KIND=JPRM),                   INTENT(IN)  :: ALPHA
!    REAL(KIND=JPRL), DIMENSION(:,:,:), INTENT(IN)  :: AARRAY
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEA
!    REAL(KIND=JPRL), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEB
!    REAL(KIND=JPRM),                   INTENT(IN)  :: BETA
!    REAL(KIND=JPRM), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEC
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT
!
!    CALL CUDA_STCGEMM_BATCHED( &
!      & TRANSA, TRANSB, &
!      & M, N, K, &
!      & ALPHA, &
!      & AARRAY, LDA, STRIDEA, &
!      & BARRAY, LDB, STRIDEB, &
!      & BETA, &
!      & CARRAY, LDC, STRIDEC, &
!      & BATCHCOUNT)
!  END SUBROUTINE CUDA_STCGEMM_BATCHED_OVERLOAD




END MODULE CUDA_GEMM_BATCHED_MOD
