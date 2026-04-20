! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE

SUBROUTINE DIR_TRANS_FIELD_API(KRESOL,               &
                             & YDFSCALAR,YDFU, YDFV, &
                             & YDFSPSCALAR, YDFSPVOR,YDFSPDIV)

!**** *DIR_TRANS_FIELD_API* - Field API interface to direct spectral transform

!     Purpose.
!     --------
!        Allow to call DIR_TRANS with a list of fields from field API

!**   Interface.
!     ----------
!     CALL DIR_TRANS_FIELD_API(...)

!     Explicit arguments :
!     --------------------
!  input
!       KRESOL         - Resolution tag
!       YDFSCALAR(:)   - List of grid-point scalar fields
!       YDFU(:)        - List of grid-point vector fields (u)
!       YDFV(:)        - List of grid-point vector fields (v)
!  output
!       YDFSPSCALAR(:) - List of spectral scalar fields
!       YDFSPVOR(:)    - List of spectral vector fields (vorticity)
!       YDFSPDIV(:)    - List of spectral vector fields (divergence)
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE ECTRANS_FIELD_API_MOD, ONLY: FIELD_GRID, FIELD_SPEC
USE PARKIND1  ,ONLY : JPIM     ,JPRB


TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)
TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)

TYPE(FIELD_SPEC),INTENT(INOUT), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)
TYPE(FIELD_SPEC),INTENT(INOUT), OPTIONAL  :: YDFSPSCALAR(:)

INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KRESOL

END SUBROUTINE DIR_TRANS_FIELD_API
END INTERFACE
