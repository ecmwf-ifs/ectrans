SUBROUTINE TRANS_INQ(KRESOL,KSPEC,KSPEC2,KSPEC2G,KSPEC2MX,KNUMP,&
                    &KGPTOT,KGPTOTG,KGPTOTMX,&
                    &PMU,PGW)

!**** *TRANS_INQ* - Extract information from the transform package

!     Purpose.
!     --------
!     Interface routine for extracting information from the T.P.

!**   Interface.
!     ----------
!     CALL TRANS_INQ(...)
!     Explicit arguments : All arguments are optional.
!     -------------------- 
!     KRESOL   - resolution tag for which info is required ,default is the
!                first defined resulution (input)
!     KSPEC    - number of complex spectral coefficients on this PE
!     KSPEC2   - 2*KSPEC 
!     KSPEC2G  - global KSPEC2
!     KSPEC2MX - maximun KSPEC2 among all PEs
!     KNUMP    - Number of spectral waves handled by this PE
!     KGPTOT   - Total number of grid columns on this PE
!     KGPTOTG  - Total number of grid columns on the Globe
!     KGPTOTMX - Maximum number of grid columns on any of the PEs
!     PMU      - sin(Gaussian latitudes)
!     PGW      - Gaussian weights

! 
!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL

INTEGER_M ,OPTIONAL, INTENT(OUT) :: KSPEC
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KSPEC2
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KSPEC2G
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KSPEC2MX
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KNUMP
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KGPTOT
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KGPTOTG
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KGPTOTMX

REAL_B    ,OPTIONAL, INTENT(OUT) :: PMU(:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PGW(:)


END SUBROUTINE TRANS_INQ






