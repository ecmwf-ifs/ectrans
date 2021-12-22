      SUBROUTINE SYMINV(A,NDIM,N,COND,V)
!
!**** *SYMINV*   - Inversion of a symmetric positive definite matrix.
!
!     Purpose.    Inversion of a symmetric positive definite matrix.
!     --------
!     Invert in place the lower triangle of A (i.e. A(i,j) i >= j)
!     A is a symmetric positive definite matrix.
!     The upper triangle of A (i.e. A(i,j) i<j) is not used nor altered.
!
!**   Interface.
!     ----------
!        *CALL* *SYMINV( ... )
!
!        Explicit arguments :
!        --------------------
!         A:       - in:  symmetric matrix to be inverted        (inout)
!                    out: its inverse
!         NDIM:    - dimension for matrices and arrays           (in)
!         N:       - actual dimension                            (in)
!         COND:    - matrix conditioning                         (out)
!         V:       - work array                                  (out)
!         DEPS:    - "epsilon" for testing positive definiteness (in-opt)
!
!        Implicit arguments :
!        --------------------
!        none.
!
!     Method.
!     -------
!        See documentation
!
!     Externals.
!     ----------
!        None.
!
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!
!     Author.
!     -------
!        ??? (former SYMINV of ECLIB).
!
!     Modifications.
!     --------------
!      K. Yessad (Aug 2009): improve comments.
!     F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------
!
      USE PARKIND1, ONLY : JPRD, JPIM, JPRB
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: NDIM
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: I,IMX,J,K,NDIA,NM1
      REAL(KIND=JPRB) :: A(NDIM,NDIM)
      REAL(KIND=JPRB) :: COND
      REAL(KIND=JPRB) :: V(NDIM)
      REAL(KIND=JPRB) :: DEPS
      REAL(KIND=JPRB) :: EPS,S,X,ZMX
      LOGICAL :: LLDOUBLE
      INTEGER(KIND=JPIM),EXTERNAL :: IDAMAX, ISAMAX
!
!     ------------------------------------------------------------------
!
      IF(N.LT.1)RETURN

      LLDOUBLE = (JPRB == JPRD)
      EPS=0.0_JPRB
      NDIA=NDIM+1
      IF (LLDOUBLE) THEN 
        IMX=IDAMAX(N,A(1,1),NDIA)
      ELSE
        IMX=ISAMAX(N,A(1,1),NDIA)
      ENDIF
      ZMX=A(IMX,IMX)*REAL(N*N,KIND=JPRB)
!     SECTION 1 & 2 COMBINED IN THIS VERSION
!**** 1. ROOT FREE CHOLESKY DECOMPOSITION  A =L D L(TRANSPOSE)
      J=1
      IF(A(1,1).LE.EPS) GO TO 91
      V(1)=1./A(1,1)
      IF(N.EQ.1) GO TO 20
      X=A(2,1)*V(1)
      A(2,2)=A(2,2)-X*A(2,1)
      A(2,1)=X
      IF(A(2,2).LE.EPS) GO TO 91
      V(2)=1./A(2,2)
      IF(N.EQ.2) GO TO 20
      DO 14 I=3,N
      DO 12 J=3,I
      S=A(I,J-1)
      DO 11 K=3,J
 11   S=S-A(I,K-2)*A(J-1,K-2)
      A(I,J-1)=S
 12   CONTINUE
      S=A(I,I)
      DO 13 J=2,I
      X=A(I,J-1)*V(J-1)
      S=S-X*A(I,J-1)
 13   A(I,J-1)=X
      A(I,I)=S
!         CHECK FOR POSITIVE-DEFINITENESS AND INVERT DIAGONAL MATRIX D.
      IF(A(I,I).LE.EPS) GO TO 91
      V(I)=1.0_JPRB/A(I,I)
 14   CONTINUE
!
!**** 2.  COPY INVERSE OF D WHICH HAS ALREADY BEEN CALCULATED.
 20   DO 21 J=1,N
 21   A(J,J)=V(J)
      IF(N.EQ.1) GO TO 50
!
!**** 3.  INVERSION OF L
 30   A(2,1)=-A(2,1)
      NM1=N-1
      IF(N.EQ.2) GO TO 40
      DO 33 I=2,NM1
      DO 32 J=2,I
      S=A(I+1,J-1)
      DO 31 K=J,I
 31   S=S+A(I+1,K)*A(K,J-1)
 32   A(I+1,J-1)=-S
 33   A(I+1,I)=-A(I+1,I)
!
!**** 4.  INV A = INV L(TRANSPOSE) * INV D * INV L
 40   DO 44 J=2,N
      S=A(J-1,J-1)
      DO 41 I=J,N
      X=A(I,I)*A(I,J-1)
      S=S+A(I,J-1)*X
 41   A(I,J-1)=X
      A(J-1,J-1)=S
      IF(J.EQ.N) GO TO 50
      DO 43 I=J,NM1
      S=A(I,J-1)
      DO 42 K=I,NM1
 42   S=S+A(K+1,I)*A(K+1,J-1)
      A(I,J-1)=S
 43   CONTINUE
 44   CONTINUE
!
 50   IF (LLDOUBLE) THEN
        IMX=IDAMAX(N,A(1,1),NDIA)
      ELSE
        IMX=ISAMAX(N,A(1,1),NDIA)
      ENDIF
      COND=1.0_JPRB/ABS(A(IMX,IMX)*ZMX)
      RETURN
!
 91   COND=-REAL(J,KIND=JPRB)
      RETURN
      ENDSUBROUTINE SYMINV
