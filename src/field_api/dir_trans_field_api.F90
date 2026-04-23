! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DIR_TRANS_FIELD_API(KRESOL,                &
                             & YDFSCALAR, YDFU, YDFV, &
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
!      output
!       YDFSPVOR(:)    - List of spectral vector fields (vorticity)
!       YDFSPDIV(:)    - List of spectral vector fields (divergence)
!       YDFSPSCALAR(:) - List of spectral scalar fields
!      input
!       KRESOL         -
!       YDFSCALAR(:)   - List of grid-point scalar fields
!       YDFU(:)        - List of grid-point vector fields (u)
!       YDFV(:)        - List of grid-point vector fields (v)
!      output
!       YDFSPSCALAR(:) - List of spectral scalar fields
!       YDFSPVOR(:)    - List of spectral vector fields (vorticity)
!       YDFSPDIV(:)    - List of spectral vector fields (divergence)
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE ECTRANS_FIELD_API_MOD, ONLY: FIELD_GRID, FIELD_SPEC, SPEC_VIEW, GRID_VIEW, LS_COUNT, LG_COUNT, LS, LG, &
                                 & IVSET_PTR, MAKE_FIELD_VIEW                                 
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW
USE PARKIND1  ,ONLY : JPIM,JPRB, JPRD

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KRESOL
TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFSCALAR(:), YDFU(:), YDFV(:)
TYPE(FIELD_SPEC),INTENT(INOUT), OPTIONAL  :: YDFSPSCALAR(:), YDFSPVOR(:), YDFSPDIV(:)

TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGSCALAR(:)
TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGU(:), YLFVGV(:)

TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVSSCALAR(:)
TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVSVOR(:), YLFVSDIV(:)

INTEGER, PARAMETER :: FIELD_TYPE = 0

INTEGER(KIND=JPIM) :: IFIELD
INTEGER(KIND=JPIM) :: NFIELD_UV, NFIELD_SCALAR
INTEGER(KIND=JPIM) :: IRESOL

#include "dir_trans_field_view.h"

REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',0,ZHOOK_HANDLE)

! -------------------------------------------------------------------------------------------------------------
! Convert FIELD_API input fields into FIELD_VIEW, an intermediate representation of fields specific to ectrans

NFIELD_UV = 0
NFIELD_SCALAR = 0
IF (PRESENT(YDFU) .AND. PRESENT(YDFV)) THEN
  NFIELD_UV = SIZE(YDFU)
  ENDIF
IF (PRESENT(YDFSCALAR)) THEN
  NFIELD_SCALAR = SIZE(YDFSCALAR)
  ENDIF
ALLOCATE(YLFVGU(NFIELD_UV))
ALLOCATE(YLFVGV(NFIELD_UV))
ALLOCATE(YLFVGSCALAR(NFIELD_SCALAR))
ALLOCATE(YLFVSVOR(NFIELD_UV))
ALLOCATE(YLFVSDIV(NFIELD_UV))
ALLOCATE(YLFVSSCALAR(NFIELD_SCALAR))

DO IFIELD=1,NFIELD_UV
  CALL MAKE_FIELD_VIEW(YLFVGU(IFIELD),   YDFU(IFIELD),     FIELD_TYPE, .FALSE., .TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVGV(IFIELD),   YDFV(IFIELD),     FIELD_TYPE, .FALSE., .TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVSVOR(IFIELD), YDFSPVOR(IFIELD), FIELD_TYPE, .FALSE., .TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVSDIV(IFIELD), YDFSPDIV(IFIELD), FIELD_TYPE, .FALSE., .TRUE.)
ENDDO
DO IFIELD=1,NFIELD_SCALAR
  CALL MAKE_FIELD_VIEW(YLFVGSCALAR(IFIELD), YDFSCALAR(IFIELD),   FIELD_TYPE, .FALSE., .FALSE.)
  CALL MAKE_FIELD_VIEW(YLFVSSCALAR(IFIELD), YDFSPSCALAR(IFIELD), FIELD_TYPE, .FALSE., .FALSE.)
ENDDO

IRESOL = 1
IF (PRESENT(KRESOL)) THEN
  IRESOL = KRESOL
ENDIF
CALL DIR_TRANS_FIELD_VIEW(IRESOL, YLFVGSCALAR, YLFVGU, YLFVGV, YLFVSSCALAR, YLFVSVOR, YLFVSDIV)

DEALLOCATE(YLFVGU)
DEALLOCATE(YLFVGV)
DEALLOCATE(YLFVGSCALAR)
DEALLOCATE(YLFVSVOR)
DEALLOCATE(YLFVSDIV)
DEALLOCATE(YLFVSSCALAR)

IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',1,ZHOOK_HANDLE)

END SUBROUTINE DIR_TRANS_FIELD_API
