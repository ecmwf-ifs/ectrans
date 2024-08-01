! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE DIST_GRID_32(PGPG,KPROMA,KFDISTG,KFROM,KRESOL,PGP)

! begin_doc_block
! ## `DIST_GRID_32`
!
! @note
! This subroutine is deprecated and will be removed in a future release.
! @endnote
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB ,JPRM


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRM)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRM)             , INTENT(OUT) :: PGP(:,:,:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID_32

END INTERFACE
