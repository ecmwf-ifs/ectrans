MODULE CPLEDN_MOD
CONTAINS
SUBROUTINE CPLEDN(KN,KDBLE,PX,PDX,KFLAG,PW,PXN,PDXN,PXMOD)

!**** *CPLEDN* - Routine to compute the Legendre polynomial of degree N

!     Purpose.
!     --------
!           Computes Legendre polynomial of degree N

!**   Interface.
!     ----------
!        *CALL* *CPLEDN(KN,KDBLE,PX,PDX,KFLAG,PW,PXN,PDXN,PXMOD)*

!        Explicit arguments :
!        --------------------
!          KN       :  Degree of the Legendre polynomial              (in)
!          KDBLE    :  0, single precision                            (in)
!                      1, double precision
!          PX       :  abcissa where the computations are performed   (in)
!          PDX      :  id in double precision                         (in)
!          KFLAG    :  When KFLAG.EQ.1 computes the weights           (in)
!          PW       :  Weight of the quadrature at PXN                (out)
!          PXN      :  new abscissa (Newton iteration)                (out)
!          PDXN     :  id in double precision                         (out)
!          PXMOD    :  PXN-PX                                         (out)

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        Michel Rochas, 90-08-30 (Lobatto+cleaning)
!        K. Yessad (Sep 2008): cleaning, improve comments.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARKIND2  ,ONLY : JPRH

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KN
INTEGER(KIND=JPIM),INTENT(IN) :: KDBLE
REAL(KIND=JPRB),INTENT(IN)    :: PX
REAL(KIND=JPRH),INTENT(IN)    :: PDX
INTEGER(KIND=JPIM),INTENT(IN) :: KFLAG
REAL(KIND=JPRB),INTENT(OUT)   :: PW
REAL(KIND=JPRB),INTENT(OUT)   :: PXN
REAL(KIND=JPRH),INTENT(OUT)   :: PDXN
REAL(KIND=JPRB),INTENT(OUT)   :: PXMOD

!     ------------------------------------------------------------------

REAL(KIND=JPRH) :: ZDLX,ZDLK,ZDLKM1,ZDLKM2,ZDLLDN,ZDLXN,ZDLMOD

INTEGER(KIND=JPIM), PARAMETER :: JPKS=KIND(PX)
INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(PDX)

INTEGER(KIND=JPIM) :: IZN, JN
REAL(KIND=JPRB) :: ZK, ZKM1, ZKM2, ZLDN, ZMOD, ZX, ZXN

!      -----------------------------------------------------------------

!*       1. NEWTON ITERATION STEP.
!           ----------------------

IZN = KN

ZK = 0.0_JPRB
ZDLK = 0.0_JPRB
ZDLXN = 0.0_JPRB

IF(KDBLE == 0)THEN

  !*       1.1   "JPRB" PRECISION.

  ZKM2 = 1
  ZKM1 = PX
  ZX   = PX
  DO JN=2,IZN
    ZK = (REAL(2*JN-1,JPRB)*ZX*ZKM1-REAL(JN-1,JPRB)*ZKM2)/REAL(JN,JPRB)
    ZKM2 = ZKM1
    ZKM1 = ZK
  ENDDO
  ZKM1 = ZKM2
  ZLDN = (REAL(KN,JPRB)*(ZKM1-ZX*ZK))/(1.0_JPRB-ZX*ZX)
  ZMOD = -ZK/ZLDN
  ZXN = ZX+ZMOD
  PXN = ZXN
  PDXN = REAL(ZXN,JPKD)
  PXMOD = ZMOD

ELSE

  !*       1.2   "JPRH" PRECISION.

  ZDLKM2 = 1.0_JPRB
  ZDLKM1 = PDX
  ZDLX = PDX
  DO JN=2,IZN
    ZDLK = (REAL(2*JN-1,JPKD)*ZDLX*ZDLKM1-REAL(JN-1,JPKD)*ZDLKM2)/REAL(JN,JPKD)
    ZDLKM2 = ZDLKM1
    ZDLKM1 = ZDLK
  ENDDO
  ZDLKM1 = ZDLKM2
  ZDLLDN = (REAL(KN,JPKD)*(ZDLKM1-ZDLX*ZDLK))/(1.0_JPRB-ZDLX*ZDLX)
  ZDLMOD = -ZDLK/ZDLLDN
  ZDLXN = ZDLX+ZDLMOD
  PXN = REAL(ZDLXN,JPKS)
  PDXN = ZDLXN
  PXMOD = REAL(ZDLMOD,JPKS)
ENDIF

!     ------------------------------------------------------------------

!*         2.    Computes weight.
!                ----------------

IF(KFLAG == 1)THEN
  ZDLKM2 = 1.0_JPRB
  ZDLKM1 = ZDLXN
  ZDLX = ZDLXN
  DO JN=2,IZN
    ZDLK = (REAL(2*JN-1,JPKD)*ZDLX*ZDLKM1-REAL(JN-1,JPKD)*ZDLKM2)/REAL(JN,JPKD)
    ZDLKM2 = ZDLKM1
    ZDLKM1 = ZDLK
  ENDDO
  ZDLKM1 = ZDLKM2
  PW = REAL((1.0_JPRB-ZDLX*ZDLX)/(REAL(KN*KN,JPKD)*ZDLKM1*ZDLKM1),JPKS)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE CPLEDN
END MODULE CPLEDN_MOD
