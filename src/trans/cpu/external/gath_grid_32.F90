! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! Single-precision version of DIST_GRID.
! This subroutine no longer has a purpose and is deprecated. It will be removed in a future release.
SUBROUTINE GATH_GRID_32(PGPG,KPROMA,KFGATHG,KTO,KRESOL,PGP)

USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
USE PARKIND1, ONLY: JPIM ,JPRM

IMPLICIT NONE

REAL(KIND=JPRM)    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRM)             , INTENT(IN)  :: PGP(:,:,:)

!     ------------------------------------------------------------------

CALL ABORT_TRANS('GATH_GRID_32: This subroutine is deprecated.')

!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID_32

