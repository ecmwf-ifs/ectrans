SUBROUTINE GATH_SPEC(PSPECG,KFDISTG,KTO,KVSET,KRESOL,PSPEC)

!**** *GATH_SPEC* - Gather global spectral array from processors

!     Purpose.
!     --------
!        Interface routine for gathering spectral array

!**   Interface.
!     ----------
!     CALL GATH_SPEC(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be gathered
!     KTO(:)      - Processor responsible for gathering each field
!     KVSET(:)    - "B-Set" for each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PSPEC(:,:)  - Local spectral array
!
!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(OUT)  :: PSPECG(:,:)
INTEGER_M          , INTENT(IN)  :: KFDISTG
INTEGER_M          , INTENT(IN)  :: KTO(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL_B    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)


!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC

