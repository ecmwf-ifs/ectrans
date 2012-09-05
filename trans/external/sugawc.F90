SUBROUTINE SUGAWC(KDGLG,PMU,PW)

!**** *SUGAWC* - Compute Gaussian latitudes and weights 

!     Purpose.
!     --------
!     Compute Gaussian latitudes and weights.

!**   Interface.
!     ----------
!     CALL SUGAWC(...)

!     Explicit arguments :
!     -------------------- 
!      INPUT:
!       KDGLG    - number of latitudes.

!      OUTPUT:
!       PMU      - sine of Gaussian latitudes.
!       PW       - Gaussian weights.

!     Method.
!     -------

!     Externals.  SUGAW
!     ----------  

!     Author.
!     -------
!        K. Yessad, from SUGAWA and SULEG (trans)
!        Original : May 2012

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE SUGAW_MOD

!endif INTERFACE

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)  :: KDGLG
REAL(KIND=JPRB)    ,INTENT(OUT) :: PMU(:)
REAL(KIND=JPRB)    ,INTENT(OUT) :: PW(:)

!ifndef INTERFACE

REAL(KIND=JPRB)    :: ZANM
INTEGER(KIND=JPIM) :: ISTART,IK,IODD,JN,JGL
REAL(KIND=JPRB) :: ZFN(0:KDGLG,0:KDGLG)
REAL(KIND=JPRB) :: ZFNN

!     ------------------------------------------------------------------

! * preliminary calculations to compute input quantities ZANM and ZFN
!   (k.y.: coded after what I found in tfl/module/suleg_mod.F90).
ISTART=1
! Belousov, Swarztrauber use ZFN(0,0)=SQRT(2._JPRB)
! IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
ZFN(0,0)=2._JPRB
DO JN=ISTART,KDGLG
  ZFNN=ZFN(0,0)
  DO JGL=1,JN
    ZFNN=ZFNN*SQRT(1._JPRB-0.25_JPRB/REAL(JGL**2,JPRB))
  ENDDO
  IODD=MOD(JN,2)
  ZFN(JN,JN)=ZFNN
  DO JGL=2,JN-IODD,2
    ZFN(JN,JN-JGL)=ZFN(JN,JN-JGL+2)*REAL((JGL-1)*(2*JN-JGL+2),JPRB)/REAL(JGL*(2*JN-JGL+1),JPRB)
  ENDDO
ENDDO

ZANM=SQRT(REAL(2*KDGLG+1,JPRB)*REAL(KDGLG**2,JPRB)/REAL(2*KDGLG-1,JPRB))

! * call to SUGAW (output: PW, PMU):
CALL SUGAW(KDGLG,0,KDGLG,PMU,PW,ZANM,ZFN)

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE SUGAWC

