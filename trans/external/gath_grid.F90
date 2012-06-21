SUBROUTINE GATH_GRID(PGPG,KPROMA,KFGATHG,KTO,KRESOL,PGP)

!**** *GATH_GRID* - Gather global gridpoint array from processors

!     Purpose.
!     --------
!        Interface routine for gathering gripoint array

!**   Interface.
!     ----------
!     CALL GATH_GRID(...)

!     Explicit arguments : 
!     -------------------- 
!     PGPG(:,:)   - Global gridpoint array
!     KFGATHG     - Global number of fields to be gathered
!     KPROMA      - blocking factor for gridpoint input
!     KTO(:)      - Processor responsible for gathering each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PGP(:,:,:)  - Local spectral array
!
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  GATH_GRID_CTL -  control routine

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
USE GATH_GRID_CTL_MOD
USE ABORT_TRANS_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KFGATHG
INTEGER_M          , INTENT(IN)  :: KTO(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B             , INTENT(IN)  :: PGP(:,:,:)

!ifndef INTERFACE

INTEGER_M :: IFRECV,J,IUBOUND(3),IPROMA,IGPBLKS

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

IPROMA = D%NGPTOT
IF(PRESENT(KPROMA)) THEN
  IPROMA = KPROMA
ENDIF
IGPBLKS = (D%NGPTOT-1)/IPROMA+1


IF(UBOUND(KTO,1) < KFGATHG) THEN
 CALL ABORT_TRANS('GATH_GRID: KTO TOO SHORT!')
ENDIF
 
IFRECV = 0
DO J=1,KFGATHG
  IF(KTO(J) < 1 .OR. KTO(J) > NPROC) THEN
    WRITE(NERR,*) 'GATH_GRID:ILLEGAL KTO VALUE',KTO(J),J
    CALL ABORT_TRANS('GATH_GRID:ILLEGAL KTO VALUE')
  ENDIF
  IF(KTO(J) == MYPROC) IFRECV = IFRECV+1
ENDDO

IUBOUND=UBOUND(PGP)
IF(IUBOUND(1) < IPROMA) THEN
  WRITE(NOUT,*)'GATH_GRID:FIRST DIM. OF PGP TOO SMALL ',IUBOUND(1),IPROMA
  CALL ABORT_TRANS('GATH_GRID:FIRST DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(2) < KFGATHG) THEN
  WRITE(NOUT,*)'GATH_GRID:SEC. DIM. OF PGP TOO SMALL ',IUBOUND(2),KFGATHG
  CALL ABORT_TRANS('GATH_GRID:SECOND DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(3) < IGPBLKS) THEN
  WRITE(NOUT,*)'GATH_GRID:THIRD DIM. OF PGP TOO SMALL ',IUBOUND(3),IGPBLKS
  CALL ABORT_TRANS('GATH_GRID:THIRD DIMENSION OF PGP TOO SMALL ')
ENDIF

IF(IFRECV > 0) THEN
  IF(.NOT.PRESENT(PGPG)) THEN
    CALL ABORT_TRANS('GATH_GRID:PGPG MISSING')
  ENDIF
  IF(UBOUND(PGPG,1) < D%NGPTOTG) THEN
    CALL ABORT_TRANS('GATH_GRID:FIRST DIMENSION OF PGPG TOO SMALL')
  ENDIF 
 IF(UBOUND(PGPG,2) < IFRECV) THEN
    CALL ABORT_TRANS('GATH_GRID:SECOND DIMENSION OF PGPG TOO SMALL')
  ENDIF
ENDIF

CALL GATH_GRID_CTL(PGPG,KFGATHG,IPROMA,KTO,PGP)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID

