SUBROUTINE SPECNORM(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)

!**** *SPECNORM* - Compute global spectral norms

!     Purpose.
!     --------
!        Interface routine for computing spectral norms

!**   Interface.
!     ----------
!     CALL SPECNORM(...)

!     Explicit arguments : All arguments optional
!     -------------------- 
!     PSPEC(:,:)  - Spectral array
!     KVSET(:)    - "B-Set" for each field
!     KMASTER     - processor to recieve norms
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PMET(:)     - metric
!     PNORM(:)    - Norms (output for processor KMASTER)
!
!     Method.
!     -------

!     Externals.  SET_RESOL - set resolution
!     ----------  SPNORM_CTL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Declaration of arguments


REAL_B    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KMASTER
REAL_B    ,OPTIONAL, INTENT(IN)  :: PMET(:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PNORM(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL

!     ------------------------------------------------------------------

END SUBROUTINE SPECNORM

