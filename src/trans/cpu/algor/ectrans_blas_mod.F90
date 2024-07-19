! (C) Copyright 2024- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!====================================================================
MODULE ECTRANS_BLAS_MOD
!====================================================================
! Author: Willem Deconinck (ECMWF)
!
! This module provides interfaces for BLAS  routines such as
! DGEMM/SGEMM and DGEMV/SGEMV
! The correct overload is used depending on the precision of the arguments
!====================================================================


USE EC_PARKIND, ONLY : JPRD, JPRM, JPIM

IMPLICIT NONE

PRIVATE

PUBLIC :: GEMM, GEMV

!---------------------------------------------------------------------

INTERFACE GEMM
! GEMM  performs one of the matrix-matrix operations
!
!    C := alpha*op( A )*op( B ) + beta*C,
!
! where  op( X ) is one of
!
!    op( X ) = X   or   op( X ) = X**T,
!
! alpha and beta are scalars, and A, B and C are matrices, with op( A )
! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

  ! SGEMM
  MODULE PROCEDURE GEMM_SP         ! Matrix arguments as array, (alpha,beta) in SP
  MODULE PROCEDURE GEMM_SP_DP      ! Matrix arguments as array, (alpha,beta) in DP
  MODULE PROCEDURE GEMM_SCAL_SP    ! Matrix arguments as scalar (address), (alpha,beta) in SP
  MODULE PROCEDURE GEMM_SCAL_SP_DP ! Matrix arguments as scalar (address), (alpha,beta) in DP

  ! DGEMM
  MODULE PROCEDURE GEMM_DP         ! Matrix arguments as array, (alpha,beta) in DP
  MODULE PROCEDURE GEMM_DP_SP      ! Matrix arguments as array, (alpha,beta) in SP
  MODULE PROCEDURE GEMM_SCAL_DP    ! Matrix arguments as scalar (address), (alpha,beta) in DP
  MODULE PROCEDURE GEMM_SCAL_DP_SP ! Matrix arguments as scalar (address), (alpha,beta) in SP
END INTERFACE

!---------------------------------------------------------------------

INTERFACE GEMV
! GEMV  performs one of the matrix-vector operations
!
!    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!
! where alpha and beta are scalars, x and y are vectors and A is an
! m by n matrix.

  ! SGEMV
  MODULE PROCEDURE GEMV_SP         ! Matrix/Vector arguments as array, (alpha,beta) in SP
  MODULE PROCEDURE GEMV_SP_DP      ! Matrix/Vector arguments as array, (alpha,beta) in DP
  MODULE PROCEDURE GEMV_SCAL_SP    ! Matrix/Vector arguments as scalar (address), (alpha,beta) in SP
  MODULE PROCEDURE GEMV_SCAL_SP_DP ! Matrix/Vector arguments as scalar (address), (alpha,beta) in DP

  ! DGEMV
  MODULE PROCEDURE GEMV_DP         ! Matrix/Vector arguments as array, (alpha,beta) in DP
  MODULE PROCEDURE GEMV_DP_SP      ! Matrix/Vector arguments as array, (alpha,beta) in SP
  MODULE PROCEDURE GEMV_SCAL_DP    ! Matrix/Vector arguments as scalar (address), (alpha,beta) in DP
  MODULE PROCEDURE GEMV_SCAL_DP_SP ! Matrix/Vector arguments as scalar (address), (alpha,beta) in SP
END INTERFACE

!---------------------------------------------------------------------

!====================================================================
CONTAINS
!====================================================================

SUBROUTINE GEMM_SCAL_DP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A, B
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: C

  CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

END SUBROUTINE GEMM_SCAL_DP

!---------------------------------------------------------------------

SUBROUTINE GEMM_SCAL_DP_SP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A, B
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: C

  CALL GEMM_SCAL_DP(TRANSA, TRANSB, M, N, K, REAL(ALPHA,JPRD), A, LDA, B, LDB, REAL(BETA,JPRD), C, LDC)

END SUBROUTINE GEMM_SCAL_DP_SP

!---------------------------------------------------------------------

SUBROUTINE GEMM_DP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A(LDA,*), B(LDB,*)
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: C(LDC,*)

  CALL GEMM_SCAL_DP(TRANSA, TRANSB, M, N, K, ALPHA, A(1,1), LDA, B(1,1), LDB, BETA, C(1,1), LDC)

END SUBROUTINE GEMM_DP

!---------------------------------------------------------------------

SUBROUTINE GEMM_DP_SP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A(LDA,*), B(LDB,*)
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: C(LDC,*)

  CALL GEMM_SCAL_DP(TRANSA, TRANSB, M, N, K, REAL(ALPHA,JPRD), A(1,1), LDA, B(1,1), LDB, REAL(BETA,JPRD), C(1,1), LDC)

END SUBROUTINE GEMM_DP_SP

!---------------------------------------------------------------------

SUBROUTINE GEMM_SCAL_SP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE, INTRINSIC :: IEEE_EXCEPTIONS, ONLY : IEEE_GET_HALTING_MODE, IEEE_SET_HALTING_MODE, IEEE_INVALID

  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A, B
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: C

#ifdef WITH_IEEE_HALT
  LOGICAL, PARAMETER :: LL_IEEE_HALT = .TRUE.
#else
  LOGICAL, PARAMETER :: LL_IEEE_HALT = .FALSE.
#endif
  LOGICAL :: LL_HALT_INVALID = .FALSE.

  IF (LL_IEEE_HALT) THEN
    CALL IEEE_GET_HALTING_MODE(IEEE_INVALID,LL_HALT_INVALID)
    IF (LL_HALT_INVALID) CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)
  ENDIF

  CALL SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

  IF (LL_IEEE_HALT .AND. LL_HALT_INVALID) CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .TRUE.)

END SUBROUTINE GEMM_SCAL_SP

!---------------------------------------------------------------------

SUBROUTINE GEMM_SCAL_SP_DP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A, B
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: C

  CALL GEMM_SCAL_SP(TRANSA, TRANSB, M, N, K, REAL(ALPHA,JPRM), A, LDA, B, LDB, REAL(BETA,JPRM), C, LDC)

END SUBROUTINE GEMM_SCAL_SP_DP

!---------------------------------------------------------------------

SUBROUTINE GEMM_SP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A(LDA,*), B(LDB,*)
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: C(LDC,*)

  CALL GEMM_SCAL_SP(TRANSA, TRANSB, M, N, K, ALPHA, A(1,1), LDA, B(1,1), LDB, BETA, C(1,1), LDC)

END SUBROUTINE GEMM_SP

!---------------------------------------------------------------------

SUBROUTINE GEMM_SP_DP(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: K, LDA, LDB, LDC, M, N
  CHARACTER           ,INTENT(IN)    :: TRANSA, TRANSB
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A(LDA,*), B(LDB,*)
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: C(LDC,*)

  CALL GEMM_SCAL_SP(TRANSA, TRANSB, M, N, K, REAL(ALPHA,JPRM), A(1,1), LDA, B(1,1), LDB, REAL(BETA,JPRM), C(1,1), LDC)

END SUBROUTINE GEMM_SP_DP

!---------------------------------------------------------------------

SUBROUTINE GEMV_SCAL_SP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A, X
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: Y

  CALL SGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)

END SUBROUTINE GEMV_SCAL_SP

!---------------------------------------------------------------------

SUBROUTINE GEMV_SCAL_SP_DP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A, X
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: Y

  CALL GEMV_SCAL_SP(TRANS, M, N, REAL(ALPHA,JPRM), A, LDA, X, INCX, REAL(BETA,JPRM), Y, INCY)

END SUBROUTINE GEMV_SCAL_SP_DP

!---------------------------------------------------------------------

SUBROUTINE GEMV_SP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A(:,:), X(:)
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: Y(:)

  CALL GEMV_SCAL_SP(TRANS, M, N, ALPHA, A(1,1), LDA, X(1), INCX, BETA, Y(1), INCY)

END SUBROUTINE GEMV_SP

!---------------------------------------------------------------------

SUBROUTINE GEMV_SP_DP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRM)     ,INTENT(IN)    :: A(:,:), X(:)
  REAL(KIND=JPRM)     ,INTENT(INOUT) :: Y(:)

  CALL GEMV_SCAL_SP(TRANS, M, N, REAL(ALPHA,JPRM), A(1,1), LDA, X(1), INCX, REAL(BETA,JPRM), Y(1), INCY)

END SUBROUTINE GEMV_SP_DP

!---------------------------------------------------------------------

SUBROUTINE GEMV_SCAL_DP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A, X
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: Y

  CALL DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)

END SUBROUTINE GEMV_SCAL_DP

!---------------------------------------------------------------------

SUBROUTINE GEMV_SCAL_DP_SP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A, X
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: Y

  CALL GEMV_SCAL_DP(TRANS, M, N, REAL(ALPHA,JPRD), A, LDA, X, INCX, REAL(BETA,JPRD), Y, INCY)

END SUBROUTINE GEMV_SCAL_DP_SP

!---------------------------------------------------------------------

SUBROUTINE GEMV_DP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRD)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A(:,:), X(:)
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: Y(:)

  CALL GEMV_SCAL_DP(TRANS, M, N, ALPHA, A(1,1), LDA, X(1), INCX, BETA, Y(1), INCY)

END SUBROUTINE GEMV_DP

!---------------------------------------------------------------------

SUBROUTINE GEMV_DP_SP(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  REAL(KIND=JPRM)     ,INTENT(IN)    :: ALPHA, BETA
  INTEGER(KIND=JPIM)  ,INTENT(IN)    :: LDA, INCX, INCY, M, N
  CHARACTER           ,INTENT(IN)    :: TRANS
  REAL(KIND=JPRD)     ,INTENT(IN)    :: A(:,:), X(:)
  REAL(KIND=JPRD)     ,INTENT(INOUT) :: Y(:)

  CALL GEMV_SCAL_DP(TRANS, M, N, REAL(ALPHA,JPRD), A(1,1), LDA, X(1), INCX, REAL(BETA,JPRD), Y(1), INCY)

END SUBROUTINE GEMV_DP_SP

!====================================================================

END MODULE ECTRANS_BLAS_MOD

