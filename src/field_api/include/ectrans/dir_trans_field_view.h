! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE

SUBROUTINE DIR_TRANS_FIELD_VIEW(    &
    & YDGPSCALAR, YDGPU, YDGPV,     &
    & YDSPSCALAR, YDSPVOR, YDSPDIV, &
    & KRESOL)

USE EC_PARKIND, ONLY : JPIM
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW

! Arguments
TYPE(FIELD_VIEW),   INTENT(IN)           :: YDGPSCALAR(:), YDGPU(:), YDGPV(:)
TYPE(FIELD_VIEW),   INTENT(INOUT)        :: YDSPSCALAR(:), YDSPVOR(:), YDSPDIV(:)
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KRESOL

END SUBROUTINE

END INTERFACE