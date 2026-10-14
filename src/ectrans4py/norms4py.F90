! Spectral and grid-point norms (distributed). Single precision. Field-first inputs
! consistent with the transform wrappers: spectral (KSPEC2, KFLD) [coeff-first in
! Fortran], grid (KGPTOT, KFLD). Norms are global (gathered across tasks).
! KRESOL, the last argument of both routines, is the resolution tag returned by
! SETUP_TRANS_4PY, and the norms are computed on that resolution. A KRESOL that is not a
! resolution set up by SETUP_TRANS_4PY returns error code 101 (RESOL_CHECK4PY_MOD).
! Grid-point arrays that ecTrans reads or writes whole are declared (KGPTOT, KFLD, 1),
! ecTrans's NPROMA-blocked layout with a single block: the same memory as (KGPTOT, KFLD),
! so they are passed to ecTrans without a copy.

SUBROUTINE SPECNORM4PY(KRETURNCODE, KSPEC2, KFLD, PSPEC, PNORM, KRESOL)
! Global spectral L2 norm per field. Cf. SPECNORM.
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM, JPRB
USE RESOL_CHECK4PY_MOD, ONLY: RESOL_STATUS
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KSPEC2, KFLD, KRESOL
REAL(KIND=JPRB), DIMENSION(KSPEC2, KFLD), INTENT(IN)  :: PSPEC
REAL(KIND=JPRB), DIMENSION(KFLD),         INTENT(OUT) :: PNORM
REAL(KIND=JPRB) :: ZSP(KFLD, KSPEC2)
INTEGER(KIND=JPIM) :: IVSET(KFLD), JL
#include "specnorm.h"
KRETURNCODE = RESOL_STATUS(KRESOL)
IF (KRETURNCODE /= 0) RETURN
IVSET(:) = 1                       ! single V-set (NPRTRV=1)
DO JL = 1, KFLD
  ZSP(JL, :) = PSPEC(:, JL)
ENDDO
PNORM(:) = 0._JPRB
CALL SPECNORM(PNORM=PNORM, PSPEC=ZSP, KVSET=IVSET, KMASTER=1, KRESOL=INT(KRESOL,JPIM))
END SUBROUTINE SPECNORM4PY


SUBROUTINE GPNORM_TRANS4PY(KRETURNCODE, KGPTOT, KFLD, PGP, PAVE, PMIN, PMAX, KRESOL)
! Global grid-point average / min / max per field. Cf. GPNORM_TRANS.
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM, JPRB
USE RESOL_CHECK4PY_MOD, ONLY: RESOL_STATUS
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KGPTOT, KFLD, KRESOL
REAL(KIND=JPRB), DIMENSION(KGPTOT, KFLD, 1), INTENT(IN)  :: PGP
REAL(KIND=JPRB), DIMENSION(KFLD),         INTENT(OUT) :: PAVE, PMIN, PMAX
#include "gpnorm_trans.h"
KRETURNCODE = RESOL_STATUS(KRESOL)
IF (KRETURNCODE /= 0) RETURN
PMIN(:) = 0._JPRB
PMAX(:) = 0._JPRB
CALL GPNORM_TRANS(PGP, INT(KFLD,JPIM), INT(KGPTOT,JPIM), PAVE, PMIN, PMAX, .FALSE., &
 & KRESOL=INT(KRESOL,JPIM))
END SUBROUTINE GPNORM_TRANS4PY
