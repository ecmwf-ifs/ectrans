SUBROUTINE DIST_GRID(PGPG,KPROMA,KFDISTG,KFROM,KRESOL,PGP)

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
!     KPROMA      - required blocking factor for gridpoint input
!     KFROM(:)    - Processor resposible for distributing each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PGP(:,:)  - Local spectral array
!
!     Method.
!     -------

!     Externals.  SET_RESOL      - set resolution
!     ----------  DIST_GRID_CTL  - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET_RESOL_MOD
USE DIST_GRID_CTL_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KFROM(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B             , INTENT(OUT) :: PGP(:,:,:)

!ifndef INTERFACE

INTEGER_M :: IFSEND,J,IUBOUND(3),IPROMA,IGPBLKS

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

IPROMA = D%NGPTOT
IF(PRESENT(KPROMA)) THEN
  IPROMA = KPROMA
ENDIF
IGPBLKS = (D%NGPTOT-1)/IPROMA+1

IF(UBOUND(KFROM,1) < KFDISTG) THEN
 CALL ABOR1('DIST_GRID: KFROM TOO SHORT!')
ENDIF
 
IFSEND = 0
DO J=1,KFDISTG
  IF(KFROM(J) < 1 .OR. KFROM(J) > NPROC) THEN
    WRITE(NERR,*) 'DIST_GRID:ILLEGAL KFROM VALUE',KFROM(J),J
    CALL ABOR1('DIST_GRID:ILLEGAL KFROM VALUE')
  ENDIF
  IF(KFROM(J) == MYPROC) IFSEND = IFSEND+1
ENDDO

IUBOUND=UBOUND(PGP)
IF(IUBOUND(1) < IPROMA) THEN
  WRITE(NOUT,*)'DIST_GRID:FIRST DIM. OF PGP TOO SMALL ',IUBOUND(1),IPROMA
  CALL ABOR1('DIST_GRID:FIRST DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(2) < KFDISTG) THEN
  WRITE(NOUT,*)'DIST_GRID:SEC. DIM. OF PGP TOO SMALL ',IUBOUND(2),KFDISTG
  CALL ABOR1('DIST_GRID:SECOND DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(3) < IGPBLKS) THEN
  WRITE(NOUT,*)'DIST_GRID:THIRD DIM. OF PGP TOO SMALL ',IUBOUND(3),IGPBLKS
  CALL ABOR1('DIST_GRID:THIRD DIMENSION OF PGP TOO SMALL ')
ENDIF

IF(IFSEND > 0) THEN
  IF(.NOT.PRESENT(PGPG)) THEN
    CALL ABOR1('DIST_GRID:PGPG MISSING')
  ENDIF
  IF(UBOUND(PGPG,1) < D%NGPTOTG) THEN
    CALL ABOR1('DIST_GRID:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF 
 IF(UBOUND(PGPG,2) < IFSEND) THEN
    CALL ABOR1('DIST_GRID:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF
ENDIF


CALL DIST_GRID_CTL(PGPG,KFDISTG,IPROMA,KFROM,PGP)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID

