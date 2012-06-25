MODULE GAWL_MOD
CONTAINS
SUBROUTINE GAWL(PL,PDL,PW,PEPS,KN,KITER,PMOD)

!**** *GAWL * - Routine to perform the Newton loop

!     Purpose.
!     --------
!           Find 0 of Legendre polynomial with Newton loop
!**   Interface.
!     ----------
!        *CALL* *GAWL(PL,PDL,PW,PEPS,KN,KITER,PMOD)

!        Explicit arguments :
!        --------------------
! PL     Gaussian latitude                         (inout)
! PDL    Gaussian latitude in double precision     (out)
! PW     Gaussian weight                           (out)
! PEPS   0 of the machine                          (in)
! KN     Truncation                                (in)
! KITER  Number of iterations                      (out)
! PMOD   Last modification                         (inout)

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        Newton Loop.

!     Externals.
!     ----------
!        CPLEDN

!     Reference.
!     ----------

!     ARPEGE Documentation vol.2, ch3.

!     Author.
!     -------
!        Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-12-18
!        K. Yessad (Sep 2008): cleaning, improve comments.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARKIND2  ,ONLY : JPRH

USE CPLEDN_MOD

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB),INTENT(INOUT)  :: PL
REAL(KIND=JPRH),INTENT(OUT)    :: PDL
REAL(KIND=JPRB),INTENT(OUT)    :: PW
REAL(KIND=JPRB),INTENT(IN)     :: PEPS
INTEGER(KIND=JPIM),INTENT(IN)  :: KN
INTEGER(KIND=JPIM),INTENT(OUT) :: KITER
REAL(KIND=JPRB),INTENT(INOUT)  :: PMOD

!     ------------------------------------------------------------------

REAL(KIND=JPRH) :: ZDLX,ZDLXN
INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(ZDLX)
INTEGER(KIND=JPIM) :: IDBLE, IFLAG, ITEMAX, JTER
REAL(KIND=JPRB) :: ZW, ZX, ZXN

!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

ITEMAX = 20
ZX = PL
ZDLX = REAL(ZX,JPKD)
IDBLE = 1
IFLAG = 0

!     ------------------------------------------------------------------

!*       2. Newton iteration.
!           -----------------

DO JTER=1,ITEMAX+1
  KITER = JTER
  CALL CPLEDN(KN,IDBLE,ZX,ZDLX,IFLAG,ZW,ZXN,ZDLXN,PMOD)
  ZX = ZXN
  ZDLX = ZDLXN

  IF(IFLAG == 1) EXIT
  IF(IDBLE == 1.AND.ABS(PMOD) <= PEPS*1000._JPRB) IFLAG = 1
  IF(ABS(PMOD) <= PEPS*1000._JPRB) IDBLE = 1

ENDDO

PL = ZXN
PDL = ZDLXN
PW = ZW

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE GAWL
END MODULE GAWL_MOD


