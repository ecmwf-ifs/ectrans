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


IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KFROM(:)
REAL_B             , INTENT(OUT) :: PGP(:,:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID

