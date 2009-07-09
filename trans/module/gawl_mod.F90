MODULE GAWL_MOD
CONTAINS
SUBROUTINE GAWL(PFN,PL,DDL,PW,PEPS,KN,KITER,PMOD)

!**** *GAWL * - Routine to perform the Newton loop

!     Purpose.
!     --------
!           Find 0 of Legendre polynomial with Newton loop
!**   Interface.
!     ----------
!        *CALL* *GAWL(PFN,PL,DDL,PW,PEPS,KN,KITER,PMOD)

!        Explicit arguments :
!        --------------------
! PFN    Fourier coefficients of series expansion for the ordinary Legendre polynomials
! PL     Gaussian latitude
! DDL    Gaussian latitude in double precision
! PW     Gaussian weight
! PEPS   0 of the machine
! KN     Truncation
! KITER  Number of iterations
! PMOD   Last modification

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
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE CPLEDN_MOD

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM) :: KITER
INTEGER(KIND=JPIM) :: KN

!     DUMMY REAL SCALARS
REAL(KIND=JPRB) :: PFN(0:KN/2)
REAL(KIND=JPRB) :: PEPS
REAL(KIND=JPRB) :: PL
REAL(KIND=JPRB) :: PMOD
REAL(KIND=JPRB) :: PW

REAL(KIND=JPRB) :: DDL,DLX,DLXN

INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(DLX)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IFLAG, ITEMAX, JTER, IODD

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZW, ZX, ZXN

!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

ITEMAX = 20
ZX = PL
DLX = REAL(ZX,JPKD)
IFLAG = 0
IODD=MOD(KN,2)

!     ------------------------------------------------------------------

!*       2. Newton iteration.
!           -----------------

DO JTER=1,ITEMAX+1
  KITER = JTER
  CALL CPLEDN(KN,IODD,PFN,ZX,DLX,IFLAG,ZW,ZXN,DLXN,PMOD)
  ZX = ZXN
  DLX = DLXN
  IF(IFLAG == 1) EXIT
  IF(ABS(PMOD) <= PEPS*1000._JPRB) IFLAG = 1
ENDDO

PL = ZXN
DDL = DLXN
PW = ZW

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE GAWL
END MODULE GAWL_MOD


