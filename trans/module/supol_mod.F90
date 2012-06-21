MODULE SUPOL_MOD
CONTAINS
SUBROUTINE SUPOL(KNSMAX,DDMU,DDPOL,DDA,DDB,DDC,DDD,DDE,DDF,DDG,DDH,DDI)

!**** *SUPOL * - Routine to compute the Legendre polynomials

!     Purpose.
!     --------
!           For a given value of mu, computes the Legendre
!           polynomials.

!**   Interface.
!     ----------
!        *CALL* *SUPOL(KNSMAX,DDMU,DDPOL,DDA,DDB,DDC,DDD,DDE
!        ,DDF,DDG,DDH,DDI)

!        Explicit arguments :
!        --------------------
!              KNSMAX   :  Truncation  (triangular)
!              DDMU     :  Abscissa at which the polynomials are computed (mu)
!              DDPOL    :  Polynomials (the first index is m and the second n)


!        Implicit arguments :   None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        K. YESSAD (MAY 1998): modification to avoid underflow.
!     ------------------------------------------------------------------

#include "tsmbkind.h"
#include "hugekind.h"

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KNSMAX

REAL_H :: DDMU,DDPOL(0:KNSMAX,0:KNSMAX)
REAL_H :: DLX,DLSITA,DL1SITA,DLKM2,DLKM1,DLK,DL1,DLS
REAL_H :: DDC(0:KNSMAX,0:KNSMAX)
REAL_H :: DDD(0:KNSMAX,0:KNSMAX)
REAL_H :: DDE(0:KNSMAX,0:KNSMAX)
REAL_H :: DDA(0:KNSMAX),DDB(0:KNSMAX),DDF(0:KNSMAX)
REAL_H :: DDG(0:KNSMAX),DDH(0:KNSMAX),DDI(0:KNSMAX)

!     LOCAL INTEGER SCALARS
INTEGER_M :: JM, JN

!     LOCAL REAL SCALARS
REAL_B :: Z

!     ------------------------------------------------------------------

!*       1. First two columns.
!           ------------------

DLX=DDMU
DLSITA=SQRT(_ONE_-DLX*DLX)

! IF WE ARE LESS THAN 1Meter FROM THE POLE,
IF(ABS(REAL(DLSITA,KIND(Z))) <= SQRT(EPSILON(Z)))THEN
  DLX=1._JPRB
  DLSITA=0._JPRB
  DL1SITA=0._JPRB
ELSE
  DL1SITA=_ONE_/DLSITA
ENDIF
DLKM2=1._JPRB
DLKM1=DLX
DDPOL(0,0)=DLKM2
DDPOL(0,1)=DLKM1*DDA(1)
DDPOL(1,1)=DLSITA*DDB(1)
DO JN=2,KNSMAX
  DLK=DDF(JN)*DLX*DLKM1-DDG(JN)*DLKM2
  DL1=DDI(JN)*(DLKM1-DLX*DLK)*DL1SITA
  DDPOL(0,JN)=DLK*DDA(JN)
  DDPOL(1,JN)=DL1*DDB(JN)
  DLKM2=DLKM1
  DLKM1=DLK
ENDDO

!     ------------------------------------------------------------------

!*       2. Diagonal (the terms 0,0 and 1,1 have already been computed)
!           -----------------------------------------------------------

DLS=DL1SITA*TINY(DLS)

!OCL SCALAR
DO JN=2,KNSMAX
  DDPOL(JN,JN)=DDPOL(JN-1,JN-1)*DLSITA*DDH(JN)
  IF ( ABS(DDPOL(JN,JN))  <  DLS ) DDPOL(JN,JN)=_ZERO_
ENDDO

!     ------------------------------------------------------------------

!*       3. General recurrence.
!           -------------------

DO JN=3,KNSMAX
!DIR$ IVDEP
!OCL NOVREC
  DO JM=2,JN-1
    DDPOL(JM,JN)=DDC(JM,JN)*DDPOL(JM-2,JN-2)&
     &-DDD(JM,JN)*DDPOL(JM-2,JN-1)*DLX &
     &+DDE(JM,JN)*DDPOL(JM  ,JN-1)*DLX
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SUPOL
END MODULE SUPOL_MOD


