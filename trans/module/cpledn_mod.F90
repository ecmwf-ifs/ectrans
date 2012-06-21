MODULE CPLEDN_MOD
CONTAINS
SUBROUTINE CPLEDN(KN,KDBLE,PX,DDX,KFLAG,PW,PXN,DDXN,PXMOD)

!**** *CPLEDN* - Routine to compute the Legendre polynomial of degree N

!     Purpose.
!     --------
!           Computes Legendre polynomial of degree N

!**   Interface.
!     ----------
!        *CALL* *CPLEDN(KN,KDBLE,PX,DDX,KFLAG,PW,PXN,DDXN,PXMOD)*

!        Explicit arguments :
!        --------------------
!              KN       :  Degree of the Legendre polynomial
!              KDBLE    :  0, single precision
!                          1, double precision
!              PX       :  abcissa where the computations are performed
!              DDX      :  id in double precision
!              KFLAG    :  When KFLAG.EQ.1 computes the weights
!              PW       :  Weight of the quadrature at PXN
!              PXN      :  new abscissa (Newton iteration)
!              DDXN     :  id in double precision
!              PXMOD    :  PXN-PX

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
!     ------------------------------------------------------------------



USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARKIND2  ,ONLY : JPRH

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM) :: KDBLE
INTEGER(KIND=JPIM) :: KFLAG
INTEGER(KIND=JPIM) :: KN

!     DUMMY REAL SCALARS
REAL(KIND=JPRB) :: PW
REAL(KIND=JPRB) :: PX
REAL(KIND=JPRB) :: PXMOD
REAL(KIND=JPRB) :: PXN


REAL(KIND=JPRH) :: DDX,DDXN,DLX,DLK,DLKM1,DLKM2,DLLDN,DLXN,DLMOD
REAL(KIND=JPRH) :: DLG,DLGDN

INTEGER(KIND=JPIM), PARAMETER :: JPKS=KIND(PX)
INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(DDX)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IZN, JN

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZG, ZGDN, ZK, ZKM1, ZKM2, ZLDN, ZMOD, ZX, ZXN


!      -----------------------------------------------------------------

!*       1. Single precision computations.
!           ------------------------------

IZN = KN

ZK = 0.0_JPRB
DLK = 0.0_JPRB
DLXN = 0.0_JPRB
IF(KDBLE == 0)THEN

!*       1.1   NEWTON ITERATION STEP.

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
  DDXN = REAL(ZXN,JPKD)
  PXMOD = ZMOD

!     ------------------------------------------------------------------

!*         2.    Double precision computations.
!                ------------------------------

ELSE

!*       2.1   NEWTON ITERATION STEP.

  DLKM2 = 1.0_JPRB
  DLKM1 = DDX
  DLX = DDX
  DO JN=2,IZN
    DLK = (REAL(2*JN-1,JPKD)*DLX*DLKM1-REAL(JN-1,JPKD)*DLKM2)/REAL(JN,JPKD)
    DLKM2 = DLKM1
    DLKM1 = DLK
  ENDDO
  DLKM1 = DLKM2
  DLLDN = (REAL(KN,JPKD)*(DLKM1-DLX*DLK))/(1.0_JPRB-DLX*DLX)
  DLMOD = -DLK/DLLDN
  DLXN = DLX+DLMOD
  PXN = REAL(DLXN,JPKS)
  DDXN = DLXN
  PXMOD = REAL(DLMOD,JPKS)
ENDIF
!     ------------------------------------------------------------------

!*         3.    Computes weight.
!                ----------------


IF(KFLAG == 1)THEN
  DLKM2 = 1.0_JPRB
  DLKM1 = DLXN
  DLX = DLXN
  DO JN=2,IZN
    DLK = (REAL(2*JN-1,JPKD)*DLX*DLKM1-REAL(JN-1,JPKD)*DLKM2)/REAL(JN,JPKD)
    DLKM2 = DLKM1
    DLKM1 = DLK
  ENDDO
  DLKM1 = DLKM2
  PW = REAL((1.0_JPRB-DLX*DLX)/(REAL(KN*KN,JPKD)*DLKM1*DLKM1),JPKS)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE CPLEDN
END MODULE CPLEDN_MOD
