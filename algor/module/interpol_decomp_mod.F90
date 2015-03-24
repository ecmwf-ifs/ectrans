MODULE INTERPOL_DECOMP_MOD

! Compute Interpolative Decomposions (ID)

! See Cheng,H., Gimbutas,Z., Martinsson,P.G. and Rokhlin,V. (2005) 
! "On the compression of low rank matrices", SIAM.J.Sci.Comput.,
!  Vol. 26, No. 4, pp1389-1404

! Also lecture notes "Mulilevel compression of Linear Operators:
! Descendents of Fast Multiple Methods and Calderon-Zygmund Theory"
! P.G.Martinsson and Mark Tygert, 2011. Chapter 7.

! Author: Mats Hamrud


USE PARKIND1, ONLY : JPRD, JPIM, JPRB, JPIB
IMPLICIT NONE
CONTAINS
!===========================================================================
SUBROUTINE COMPUTE_ID(PEPS,KM,KN,PMAT,KRANK,KBCLIST,PNONIM)
IMPLICIT NONE

! Compute ID

REAL(KIND=JPRB),INTENT(IN)     :: PEPS  ! Precision for computation
                                        ! of numerical rank
INTEGER(KIND=JPIM),INTENT(IN)  :: KM    ! Number of rows in matrix pmat
INTEGER(KIND=JPIM),INTENT(IN)  :: KN    ! Number of columns in matrix pmat
REAL(KIND=JPRB)   ,INTENT(IN)  :: PMAT(:,:)  ! Original matrix
INTEGER(KIND=JPIM),INTENT(OUT) :: KRANK      ! Numerical rank
INTEGER(KIND=JPIM),INTENT(OUT) :: KBCLIST(:) ! List of  columns
REAL(KIND=JPRB)   ,INTENT(OUT) :: PNONIM(:,:)  ! Non-identity part of projection
                                               ! matrix

INTEGER(KIND=JPIM) :: JM,JN
REAL(KIND=JPRB) :: ZR(KM,KN)
REAL(KIND=JPRB),ALLOCATABLE :: ZS(:,:),ZT(:,:)
!----------------------------------------------------------------------------
!Avoid destroying input matrix
ZR(:,:) = PMAT(1:KM,1:KN)
! Householder QR
CALL ALG541(PEPS,KM,KN,ZR,KRANK,KBCLIST)

DO JN=1,KN
  DO JM=JN+1,KM
    ZR(JM,JN) = 0.0_JPRB
  ENDDO
ENDDO
! S leftmost kxk block of R
ALLOCATE(ZS(KRANK,KRANK))
DO JN=1,KRANK
  DO JM=1,KRANK
    IF(JM <= KM ) THEN
      ZS(JM,JN) = ZR(JM,JN) 
    ELSE
      ZS(JM,JN) = 0.0_JPRB
    ENDIF
  ENDDO
ENDDO
! T Rightmost kx(k-n) block of R 
ALLOCATE(ZT(KRANK,KN-KRANK))
DO JN=1,KN-KRANK
  DO JM=1,KRANK
    IF(JM <= KM ) THEN
      ZT(JM,JN) = ZR(JM,JN+KRANK) 
    ELSE
      ZT(JM,JN) = 0.0_JPRB
    ENDIF
  ENDDO
ENDDO
!Solve linear equation (BLAS level 3 routine)
IF (JPRB == JPRD) THEN
  CALL DTRSM('Left','Upper','No transpose','Non-unit',KRANK,KN-KRANK,1.0_JPRB, &
 & ZS,KRANK,ZT,KRANK)
ELSE
  CALL STRSM('Left','Upper','No transpose','Non-unit',KRANK,KN-KRANK,1.0_JPRB, &
 & ZS,KRANK,ZT,KRANK)
ENDIF

DO JM=1,KRANK
  DO JN=1,KN-KRANK
    PNONIM(JM,JN) = ZT(JM,JN)
  ENDDO
ENDDO
DEALLOCATE(ZS,ZT)
!!$IF(KRANK < KN) THEN
!!$  PRINT *,'MAXVAL PNONIM ',KM,KM,KRANK,MAXVAL( PNONIM(1:KRANK,1:KN-KRANK))
!!$ENDIF
END SUBROUTINE COMPUTE_ID

!==============================================================================
SUBROUTINE ALG541(PEPS,KM,KN,PA,KRANK,KLIST)
IMPLICIT NONE

! Householder QR with Column Pivoting
! Algorithm 5.4.1 from Matrix Computations, G.H.Golub & C.F van Loen, third ed.

! Algorithm modified to terminate at numerical precision "peps"

REAL(KIND=JPRB),INTENT(IN)     :: PEPS     ! Precision
INTEGER(KIND=JPIM),INTENT(IN)  :: KM       ! Number of rows in matrix pa
INTEGER(KIND=JPIM),INTENT(IN)  :: KN       ! Number of columns in matrix pa
REAL(KIND=JPRB),INTENT(INOUT)  :: PA(:,:)  ! On input : original matrix
                                           ! on output : R in upper triangle etc
                                           ! see Golub&Van Loen
INTEGER(KIND=JPIM),INTENT(OUT) :: KRANK    ! Numerical rank of matrix
INTEGER(KIND=JPIM),INTENT(OUT) :: KLIST(:) ! List of columns (pivots)

INTEGER(KIND=JPIM)           :: JM,JN,ISWAP,JK,IK,IIK,IM,IN,ILIST(KN)
REAL(KIND=JPRB) :: ZC(KN),ZTAU,ZSWAPA(KM),ZSWAP,ZV(KM),ZBETA,ZWORK(KN),ZTAU_IN
REAL(KIND=JPRB) :: ZTAU_REC,ZEPS
!-------------------------------------------------------------------------------
ZEPS = 10000.0_JPRB*EPSILON(ZEPS)
! Compute initial column norms,its max and the first column where c=tau
ZTAU = 0.0_JPRB
IK = 0
DO JN=1,KN
  ZC(JN) = 0.0_JPRB
  DO JM=1,KM
    ZC(JN)=ZC(JN)+PA(JM,JN)**2
  ENDDO
  IF(ZC(JN) > ZTAU) THEN
    IK = JN
    ZTAU = ZC(JN)
  ENDIF
ENDDO
ZTAU_IN = ZTAU
ZTAU_REC= ZTAU
KRANK = 0
DO WHILE (ZTAU > PEPS**2*ZTAU_IN)
  KRANK = KRANK+1
  ILIST(KRANK) = IK
! Column swap
  ZSWAPA(:) = PA(:,KRANK)
  PA(:,KRANK) = PA(:,IK)
  PA(:,IK) = ZSWAPA(:)
  ZSWAP = ZC(KRANK)
  ZC(KRANK) = ZC(IK)
  ZC(IK) = ZSWAP
! Compute Householder vector
  CALL ALG511(PA(KRANK:KM,KRANK),ZV,ZBETA)

! Apply Householder matrix
  IM = KM-KRANK+1
  IN = KN-KRANK+1
! LAPACK
  CALL DLARF('Left',IM,IN,ZV,1,ZBETA,PA(KRANK,KRANK),KM,ZWORK)
  PA(KRANK+1:KM,KRANK) = ZV(2:KM-KRANK+1)

! Update column norms
  ZTAU = 0.0_JPRB
  IF(KRANK < MIN(KM,KN)) THEN
    DO JN=KRANK+1,KN
      ZC(JN) = ZC(JN)-PA(KRANK,JN)**2
      IF(ZC(JN) > ZTAU) THEN
        IK = JN
        ZTAU = ZC(JN)
      ENDIF
    ENDDO
! Re-compute column norms when ztau < zeps*ztau_rec
!   We have disabled this test so as to recompute column nodes unconditionally.
!   This was done to resolve a SIG FPE in dtrsm @ T3999 on the XC-30.
!   IF(ZTAU < ZEPS*ZTAU_REC) THEN
      DO JN=KRANK+1,KN
        ZC(JN) = DOT_PRODUCT(PA(KRANK+1:,JN),PA(KRANK+1:,JN))
        IF(ZC(JN) > ZTAU) THEN
          IK = JN
          ZTAU = ZC(JN)
        ENDIF
      ENDDO
!!$      PRINT *,'RECOMPUTE TAU ',KRANK,ZTAU_REC,ZTAU
      ZTAU_REC = ZTAU
!   ENDIF
  ENDIF
ENDDO
! Make sure klist is filled also beyond krank
DO JN=1,KN
  KLIST(JN) = JN
ENDDO
DO JN=1,KRANK
  ISWAP = KLIST(JN)
  KLIST(JN) = KLIST(ILIST(JN))
  KLIST(ILIST(JN)) = ISWAP
ENDDO
  
END SUBROUTINE ALG541
!==============================================================================
SUBROUTINE ALG511(PX,PV,PBETA)
IMPLICIT NONE
! Compute Householder vector
! Algorithm 5.1.1 from Matrix Computations, G.H.Golub & C.F van Loen, third ed.
REAL(KIND=JPRB),INTENT(IN)     :: PX(:)
REAL(KIND=JPRB),INTENT(OUT)    :: PV(:)
REAL(KIND=JPRB),INTENT(OUT)    :: PBETA

INTEGER(KIND=JPIM) :: IN
REAL(KIND=JPRB) :: ZSIGMA,ZMU
!-------------------------------------------------------------------------------
IN = SIZE(PX)
ZSIGMA = DOT_PRODUCT(PX(2:IN),PX(2:IN))
PV(1) = 1.0_JPRB
PV(2:IN) = PX(2:IN)
IF(ZSIGMA == 0.0_JPRB) THEN
  PBETA = 0.0_JPRB
ELSE
  ZMU = SQRT(PX(1)**2+ZSIGMA)
  IF(PX(1) <= 0.0_JPRB) THEN
    PV(1) = PX(1)-ZMU
  ELSE
    PV(1) = -ZSIGMA/(PX(1)+ZMU)
  ENDIF
  PBETA = 2.0_JPRB*PV(1)**2/(ZSIGMA+PV(1)**2)
  PV(:) = PV(:)/PV(1)
ENDIF
END SUBROUTINE ALG511
!================================================================================
      
END MODULE INTERPOL_DECOMP_MOD
