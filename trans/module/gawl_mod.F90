MODULE GAWL_MOD
CONTAINS
SUBROUTINE GAWL(PL,DDL,PW,PEPS,KN,KITER,PMOD)

!**** *GAWL * - Routine to perform the Newton loop

!     Purpose.
!     --------
!           Find 0 of Legendre polynomial with Newton loop
!**   Interface.
!     ----------
!        *CALL* *GAWL(PL,DDL,PW,PEPS,KN,KITER,PMOD)

!        Explicit arguments :
!        --------------------
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
!     ------------------------------------------------------------------

#include "tsmbkind.h"
#include "hugekind.h"

USE CPLEDN_MOD

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KITER
INTEGER_M :: KN

!     DUMMY REAL SCALARS
REAL_B :: PEPS
REAL_B :: PL
REAL_B :: PMOD
REAL_B :: PW

REAL_H :: DDL,DLX,DLXN

INTEGER_M, PARAMETER :: JPKD=KIND(DLX)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IDBLE, IFLAG, ITEMAX, JTER

!     LOCAL REAL SCALARS
REAL_B :: ZW, ZX, ZXN

!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

ITEMAX = 20
ZX = PL
DLX = REAL(ZX,JPKD)
IDBLE = 1
IFLAG = 0

!     ------------------------------------------------------------------

!*       2. Newton iteration.
!           -----------------

DO JTER=1,ITEMAX+1
  KITER = JTER
  CALL CPLEDN(KN,IDBLE,ZX,DLX,IFLAG,ZW,ZXN,DLXN,PMOD)
  ZX = ZXN
  DLX = DLXN

  IF(IFLAG == 1) EXIT
  IF(IDBLE == 1.AND.ABS(PMOD) <= PEPS*1000._JPRB) IFLAG = 1
  IF(ABS(PMOD) <= PEPS*1000._JPRB) IDBLE = 1

ENDDO

PL = ZXN
DDL = DLXN
PW = ZW

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE GAWL
END MODULE GAWL_MOD


