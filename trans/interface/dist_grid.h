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
!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KFROM(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B             , INTENT(OUT) :: PGP(:,:,:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID

