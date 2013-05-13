INTERFACE
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

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)  :: KDGLG
REAL(KIND=JPRB)    ,INTENT(OUT) :: PMU(:)
REAL(KIND=JPRB)    ,INTENT(OUT) :: PW(:)

END SUBROUTINE SUGAWC

END INTERFACE
