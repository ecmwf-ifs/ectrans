SUBROUTINE DIST_SPEC(PSPECG,KFDISTG,KFROM,KVSET,KRESOL,PSPEC)

!**** *DIST_SPEC* - Distribute global spectral array among processors

!     Purpose.
!     --------
!        Interface routine for distributing spectral array

!**   Interface.
!     ----------
!     CALL DIST__SPEC(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KFROM(:)    - Processor resposible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PSPEC(:,:)  - Local spectral array
!
!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN)  :: PSPECG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KFROM(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPEC(:,:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC

