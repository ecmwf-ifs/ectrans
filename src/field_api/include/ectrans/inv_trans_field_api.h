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


SUBROUTINE INV_TRANS_FIELD_API(KRESOL,                                       &
                             & YDFSPSCALAR, YDFSPVOR,YDFSPDIV,               &
                             & YDFSCALAR, YDFU, YDFV,                        &
                             & YDFVOR,YDFDIV,                                &
                             & YDFSCALAR_NS, YDFSCALAR_EW, YDFU_EW, YDFV_EW, &
                             & KGPTOT,                                       &
                             & FSPGL_PROC)

!**** *INV_TRANS_FIELD_API* - Field API interface to inverse spectral transform

!     Purpose.
!     --------
!        Allow to call INV_TRANS with a list of fields from field API

!**   Interface.
!     ----------
!     CALL INV_TRANS_FIELD_API(...)

!     Explicit arguments :
!     --------------------
!      input
!       KRESOL           The resolution identifier
!       YDFSPSCALAR(:) - List of spectral scalar fields
!       YDFSPVOR(:)    - List of spectral vector fields (vorticity)
!       YDFSPDIV(:)    - List of spectral vector fields (divergence)
!       KGPTOT         - Number of total grid points
!       FSPGL_PROC     - procedure to be executed in fourier space
!                        before transposition

!      output
!       YDFSCALAR(:)   - List of grid-point scalar fields
!       YDFU(:)        - List of grid-point vector fields (u)
!       YDFV(:)        - List of grid-point vector fields (v)
!       YDFVOR(:)      - List of grid-point vector fields (vorticity)
!       YDFDIV(:)      - List of grid-point vector fields (divergence)
!       YDFSCALAR_NS(:) - List of grid-point scalar fields derivatives N-S
!       YDFSCALAR_EW(:) - List of grid-point scalar fields derivatives E-W
!       YDFU_EW(:)      - List of grid-point vector fields derivatives E-W (u)
!       YDFV_EW(:)      - List of grid-point vector fields derivatives E-W (v)

USE YOMHOOK, ONLY : LHOOK,   DR_HOOK, JPHOOK
USE ECTRANS_FIELD_API_MOD, ONLY: FIELD_GRID, FIELD_SPEC
USE PARKIND1, ONLY : JPIM, JPRB

#include "fspgl_intf.h"

INTEGER(KIND=JPIM),   INTENT(IN), OPTIONAL  :: KRESOL
TYPE(FIELD_SPEC),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)        ! SPECTRAL VECTOR FIELDS : VORTICITY AND DIVERGENCE FIELDS (IN)
TYPE(FIELD_SPEC),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)                  ! SPECTRAL SCALAR FIELDS (IN)

TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)                 ! GRID VECTOR FIELDS     (OUT)
TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFVOR(:),YDFDIV(:)             ! GRID VECTOR FIELDS :VORTICITY AND DIVERGENCE     (OUT)
TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)                    ! GRID SCALAR FIELDS     (OUT)

TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFU_EW(:),YDFV_EW(:)             ! GRID VECTOR FIELDS DERIVATIVES EW (OUT)
TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFSCALAR_NS(:), YDFSCALAR_EW(:)  ! GRID SCALAR FIELDS DERIVATIVES EW AND NS (OUT)

INTEGER(KIND=JPIM),   INTENT(IN)            :: KGPTOT
PROCEDURE(FSPGL_INTF), POINTER, INTENT(IN), OPTIONAL  :: FSPGL_PROC

END SUBROUTINE INV_TRANS_FIELD_API

END INTERFACE
