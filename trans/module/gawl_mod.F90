MODULE GAWL_MOD
CONTAINS
SUBROUTINE GAWL(PFN,PL,PDL,PW,PEPS,KN,KITER,PMOD)

!**** *GAWL * - Routine to perform the Newton loop

!     Purpose.
!     --------
!           Find 0 of Legendre polynomial with Newton loop
!**   Interface.
!     ----------
!        *CALL* *GAWL(PFN,PL,PDL,PW,PEPS,KN,KITER,PMOD)

!        Explicit arguments :
!        --------------------
! PFN    Fourier coefficients of series expansion
!        for the ordinary Legendre polynomials     (in)
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
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE CPLEDN_MOD

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB),INTENT(IN)     :: PFN(0:KN/2)
REAL(KIND=JPRB),INTENT(INOUT)  :: PL
REAL(KIND=JPRB),INTENT(OUT)    :: PDL
REAL(KIND=JPRB),INTENT(OUT)    :: PW
REAL(KIND=JPRB),INTENT(IN)     :: PEPS
INTEGER(KIND=JPIM),INTENT(IN)  :: KN
INTEGER(KIND=JPIM),INTENT(OUT) :: KITER
REAL(KIND=JPRB),INTENT(INOUT)  :: PMOD

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDLX,ZDLXN
INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(ZDLX)
INTEGER(KIND=JPIM) :: IFLAG, ITEMAX, JTER, IODD
REAL(KIND=JPRB) :: ZW, ZX, ZXN

!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

ITEMAX = 20
ZX = PL
ZDLX = REAL(ZX,JPKD)
IFLAG = 0
IODD=MOD(KN,2)

!     ------------------------------------------------------------------

!*       2. Newton iteration.
!           -----------------

DO JTER=1,ITEMAX+1
  KITER = JTER
  CALL CPLEDN(KN,IODD,PFN,ZX,ZDLX,IFLAG,ZW,ZXN,ZDLXN,PMOD)
  ZX = ZXN
  ZDLX = ZDLXN

  IF(IFLAG == 1) EXIT
  IF(ABS(PMOD) <= PEPS*1000._JPRB) IFLAG = 1
ENDDO

PL = ZXN
PDL = ZDLXN
PW = ZW

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE GAWL
END MODULE GAWL_MOD


