SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KLOEN,LDSPLIT,KAPSETS)

!**** *SETUP_TRANS* - Setup transform package for specific resolution

!     Purpose.
!     --------
!     To setup for making spectral transforms

!**   Interface.
!     ----------
!     CALL SETUP_TRANS(...)

!     Explicit arguments : KLOEN,LDSPLIT,KAPSETS are optional arguments
!     -------------------- 
!     KSMAX - spectral truncation required
!     KDGL  - number of gaussian latitudes
!     KLOEN(:) - number of points on each Gaussian latitude
!     LDSPLIT - true if split latitudes in grid-point space
!     KAPSETS - Number of apple sets in the distribution
!
!     KSMAX,KDGL and KLOEN are GLOBAL varaibles desribing the resolition
!     in spectral and grid-point space

!     LDSPLIT and KAPSETS describe the distribution among processors of
!     grid-point data and has no relevance if you are using a single processor
 
!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Dummy arguments

INTEGER_M ,INTENT(IN) :: KSMAX,KDGL
INTEGER_M ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER_M ,OPTIONAL,INTENT(IN) :: KAPSETS


END SUBROUTINE SETUP_TRANS


