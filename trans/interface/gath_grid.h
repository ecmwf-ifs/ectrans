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
!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KFGATHG
INTEGER_M          , INTENT(IN)  :: KTO(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B             , INTENT(IN)  :: PGP(:,:,:)


!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID

