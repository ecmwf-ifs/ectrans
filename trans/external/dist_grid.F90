SUBROUTINE DIST_GRID(PGPG,KFDISTG,KFROM,PGP)

!**** *DIST_GRID* - Distribute global gridpoint array among processors

!     Purpose.
!     --------
!        Interface routine for distributing gridpoint array

!**   Interface.
!     ----------
!     CALL DIST_GRID(...)

!     Explicit arguments : 
!     -------------------- 
!     PGPG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KFROM(:)    - Processor resposible for distributing each field
!     PGP(:,:)  - Local spectral array
!
!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE DIST_GRID_CONTROL_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KFROM(:)
REAL_B             , INTENT(OUT) :: PGP(:,:)

!ifndef INTERFACE

INTEGER_M :: IFSEND,J

!     ------------------------------------------------------------------

IF(UBOUND(KFROM,1) < KFDISTG) THEN
 CALL ABOR1('DIST_GRID: KFROM TOO SHORT!')
ENDIF
 
IFSEND = 0
DO J=1,KFDISTG
  IF(KFROM(J) == MYPROC) IFSEND = IFSEND+1
ENDDO

IF(IFSEND > 0) THEN
  IF(.NOT.PRESENT(PGPG)) THEN
    CALL ABOR1('DIST_GRID:PGPG MISSING')
  ENDIF
  IF(UBOUND(PGPG,1) < IFSEND) THEN
    CALL ABOR1('DIST_GRID:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF 
 IF(UBOUND(PGPG,2) < R%NSPEC2_G) THEN
    CALL ABOR1('DIST_GRID:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF
ENDIF


IF(UBOUND(PGP,1) < KFDISTG) THEN
  CALL ABOR1('DIST_GRID: FIRST DIMENSION OF PGP TOO SMALL')
ENDIF
IF(UBOUND(PGP,2) < D%NSPEC2 ) THEN
  CALL ABOR1('DIST_GRID: SECOND DIMENSION OF PGP TOO SMALL')
ENDIF

CALL DIST_GRID_CONTROL(PGPG,KFDISTG,KFROM,PGP)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID

