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



#include "tsmbkind.h"
#include "hugekind.h"

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KDBLE
INTEGER_M :: KFLAG
INTEGER_M :: KN

!     DUMMY REAL SCALARS
REAL_B :: PW
REAL_B :: PX
REAL_B :: PXMOD
REAL_B :: PXN


REAL_H :: DDX,DDXN,DLX,DLK,DLKM1,DLKM2,DLLDN,DLXN,DLMOD
REAL_H :: DLG,DLGDN

INTEGER_M, PARAMETER :: JPKS=KIND(PX)
INTEGER_M, PARAMETER :: JPKD=KIND(DDX)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IZN, JN

!     LOCAL REAL SCALARS
REAL_B :: ZG, ZGDN, ZK, ZKM1, ZKM2, ZLDN, ZMOD, ZX, ZXN


!      -----------------------------------------------------------------

!*       1. Single precision computations.
!           ------------------------------

IZN = KN

ZK = _ZERO_
DLK = _ZERO_
DLXN = _ZERO_
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
  ZLDN = (REAL(KN,JPRB)*(ZKM1-ZX*ZK))/(_ONE_-ZX*ZX)
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

  DLKM2 = _ONE_
  DLKM1 = DDX
  DLX = DDX
  DO JN=2,IZN
    DLK = (REAL(2*JN-1,JPKD)*DLX*DLKM1-REAL(JN-1,JPKD)*DLKM2)/REAL(JN,JPKD)
    DLKM2 = DLKM1
    DLKM1 = DLK
  ENDDO
  DLKM1 = DLKM2
  DLLDN = (REAL(KN,JPKD)*(DLKM1-DLX*DLK))/(_ONE_-DLX*DLX)
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
  DLKM2 = _ONE_
  DLKM1 = DLXN
  DLX = DLXN
  DO JN=2,IZN
    DLK = (REAL(2*JN-1,JPKD)*DLX*DLKM1-REAL(JN-1,JPKD)*DLKM2)/REAL(JN,JPKD)
    DLKM2 = DLKM1
    DLKM1 = DLK
  ENDDO
  DLKM1 = DLKM2
  PW = REAL((_ONE_-DLX*DLX)/(REAL(KN*KN,JPKD)*DLKM1*DLKM1),JPKS)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE CPLEDN
END MODULE CPLEDN_MOD
