SUBROUTINE TRANS_INQ(KRESOL,KSPEC,KSPEC2,KSPEC2MX,KNUMP,&
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

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR
USE TPM_GEOMETRY
USE TPM_FIELDS

USE SET_RESOL_MOD

!endif INTERFACE

IMPLICIT NONE

INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL

INTEGER_M ,OPTIONAL, INTENT(OUT) :: KSPEC
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KSPEC2
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KSPEC2MX
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KNUMP
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KGPTOT
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KGPTOTG
INTEGER_M ,OPTIONAL, INTENT(OUT) :: KGPTOTMX

REAL_B    ,OPTIONAL, INTENT(OUT) :: PMU(:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PGW(:)

!ifndef INTERFACE

!     ------------------------------------------------------------------


! Set current resolution
CALL SET_RESOL(KRESOL)

IF(PRESENT(KSPEC))    KSPEC    = D%NSPEC
IF(PRESENT(KSPEC2))   KSPEC2   = D%NSPEC2
IF(PRESENT(KSPEC2MX)) KSPEC2MX = D%NSPEC2MX
IF(PRESENT(KNUMP))    KNUMP    = D%NUMP
IF(PRESENT(KGPTOT))   KGPTOT   = D%NGPTOT
IF(PRESENT(KGPTOTG))  KGPTOTG  = D%NGPTOTG
IF(PRESENT(KGPTOTMX)) KGPTOTMX = D%NGPTOTMX

IF(PRESENT(PMU)) THEN
  IF(UBOUND(PMU,1) < R%NDGL) THEN
    CALL ABOR1('TRANS_INQ: PMU TOO SHORT')
  ELSE
    PMU(1:R%NDGL) = F%RMU
  ENDIF
ENDIF

IF(PRESENT(PGW)) THEN
  IF(UBOUND(PGW,1) < R%NDGL) THEN
    CALL ABOR1('TRANS_INQ: PGW TOO SHORT')
  ELSE
    PGW(1:R%NDGL) = F%RW
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE TRANS_INQ






