MODULE SUGAW_MOD
CONTAINS
SUBROUTINE SUGAW(KN,PL,DDL,PW)

#include "tsmbkind.h"
#include "hugekind.h"

USE TPM_GEN
USE GAWL_MOD

#ifdef DOC

!**** *SUGAW * - Routine to initialize the Gaussian 
!                 abcissa and the associated weights

!     Purpose.
!     --------
!           Initialize arrays PL,DDL and PW (quadrature abscissas and weights)
!**   Interface.
!     ----------
!        *CALL* *SUGAW(KN,PL,DDL,PW) *

!        Explicit arguments :
!        --------------------
!           INPUT:
!              KN       :  Number of Gauss  abscissas 
!           OUTPUT:
!              PL (KN)  :  abscissas of Gauss
!              DDL(KN)  :  idem in double precision
!              PW (KN)  :  Weights of the Gaussian integration

!     PL (i) is the abscissa i starting from the northern pole, it is
! the cosine of the colatitude of the corresponding row of the collocation
! grid. 

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!      Calls:
!      Called by SULEG.

!     Reference.
!     ----------

!     ARPEGE Documentation vol.2, ch3.

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original          : 87-10-15
!        Michel Rochas     : 90-08-30
!        Philippe Courtier : 92-12-19 Multitasking
!        Ryad El Khatib    : 94-04-20 Remove unused comdecks pardim and yomdim
!        Mats Hamrud       : 94-08-12 Printing level
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER_M,INTENT(IN) :: KN

REAL_B,INTENT(OUT) :: PL(:),PW(:)
REAL_H,INTENT(OUT) :: DDL(:)

! LOCAL REALS

REAL_B :: ZLI(KN),ZT(KN)
REAL_B :: ZREG(KN),ZMOD(KN),ZM(KN),ZR(KN)
INTEGER_M :: ITER(KN)
REAL_B :: ZD(KN),ZE(KN),ZZ(KN,KN)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IALLOW, INS2, ISYM, JGL, JROC

!     LOCAL REAL SCALARS
REAL_B :: Z, ZEPS, ZPI

LOGICAL LLP1,LLP2


!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

LLP1=.FALSE.
LLP2=.FALSE.
ZPI = _TWO_*ASIN(_ONE_)
INS2 = KN/2+MOD(KN,2)

!*       1.1 Find first approximation of the roots of the
!           Legendre polynomial of degree KN.
DO JGL=1,INS2
  Z = REAL(4*JGL-1,JPRB)*ZPI/REAL(4*KN+2,JPRB)
  PL(JGL) = COS(Z+_ONE_/(TAN(Z)*REAL(8*KN**2,JPRB)))
  ZREG(JGL) = COS(Z)
  ZLI(JGL) = PL(JGL)
ENDDO

!     ------------------------------------------------------------------

!*      2. Computes roots and weights.
!          ---------------------------

ZEPS = EPSILON(Z)
!$OMP PARALLEL DO PRIVATE(JGL)
DO JGL=1,INS2
  CALL GAWL(PL(JGL),DDL(JGL),PW(JGL),ZEPS,KN,ITER(JGL),ZMOD(JGL))
ENDDO
!$OMP END PARALLEL DO

!DIR$ IVDEP
!OCL NOVREC
DO JGL=1,KN/2
  ISYM = KN-JGL+1
  PL(ISYM) = -PL(JGL)
  DDL(ISYM) = -DDL(JGL)
  PW(ISYM) = PW(JGL)
ENDDO

!     ------------------------------------------------------------------

!*      3. Diagnostics.
!          ------------

IF(LLP2)THEN
  DO JGL=1,INS2
    ZM(JGL) = (ACOS(PL(JGL))-ACOS(ZLI(JGL)))*6371229._JPRB
    ZR(JGL) = (ACOS(PL(JGL))-ACOS(ZREG(JGL)))*6371229._JPRB
    ZT(JGL) = ACOS(PL(JGL))*180._JPRB/ZPI
  ENDDO
ENDIF

IALLOW = 10
DO JGL=1,INS2
  IF(ITER(JGL) > IALLOW)THEN
    WRITE(UNIT=NOUT,FMT='('' CONVERGENCE FAILED IN SUGAW '')')
    WRITE(UNIT=NOUT,FMT='('' ALLOWED : '',I4,''&
     &NECESSARY : '',&
     &I4)')IALLOW,ITER(JGL)
    CALL ABOR1(' FAILURE IN SUGAW ')
  ENDIF

  IF(LLP2)THEN
    WRITE(UNIT=NOUT,FMT=&
     &'('' ROW ='',I4,'' ITERATIONS='',I4,'' ROOT='',F30.20,&
     &'' WEIGHT='',F30.20,'' MODIF :'',E8.2)')JGL,ITER(JGL),PL(JGL)&
     &,PW(JGL),PL(JGL)-ZLI(JGL)
    WRITE(UNIT=NOUT,FMT=&
     &'(10X,'' LAST INC. : '',E8.2,'' MODIF IN M : '',F10.3,&
     &'' FROM THE REGULAR GRID : '',F10.3,'' COLAT '',F10.3)')&
     &ZMOD(JGL),ZM(JGL),ZR(JGL),ZT(JGL)
  ENDIF
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SUGAW
END MODULE SUGAW_MOD
