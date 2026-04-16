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
SUBROUTINE INV_TRANS_FIELD_VIEW(KRESOL, &
                               & YDSPSCALAR, YDSPVOR, YDSPDIV, &
                               & YDGPSCALAR, YDGPU, YDGPV,     &
                               & YDGPVOR,YDGPDIV, &
                               & YDGPSCALAR_NS, YDGPSCALAR_EW, YDGPU_EW, YDGPV_EW&
                               & FSPGL_PROC)

!**** *INV_TRANS_FIELD_VIEW* - Field API interface to inverse spectral transform

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
!       YDSPSCALAR(:) - List of spectral scalar fields
!       YDSPVOR(:)    - List of spectral vector fields (vorticity)
!       YDSPDIV(:)    - List of spectral vector fields (divergence)
!       FSPGL_PROC     - procedure to be executed in fourier space
!                        before transposition

!      output
!       YDGPSCALAR(:)   - List of grid-point scalar fields
!       YDGPU(:)        - List of grid-point vector fields (u)
!       YDGPV(:)        - List of grid-point vector fields (v)
!       YDGPVOR(:)      - List of grid-point vector fields (vorticity)
!       YDGPDIV(:)      - List of grid-point vector fields (divergence)
!       YDGPSCALAR_NS(:) - List of grid-point scalar fields derivatives N-S
!       YDGPSCALAR_EW(:) - List of grid-point scalar fields derivatives E-W
!       YDGPU_EW(:)      - List of grid-point vector fields derivatives E-W (u)
!       YDGPV_EW(:)      - List of grid-point vector fields derivatives E-W (v)

USE YOMHOOK, ONLY : LHOOK,   DR_HOOK, JPHOOK
USE ECTRANS_FIELD_API_BASIC_TYPE_MOD, ONLY: FIELD_BASIC_PTR
USE PARKIND1, ONLY : JPIM, JPRB

#include "fspgl_intf.h"

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)        ! SPECTRAL VECTOR FIELDS : VORTICITY AND DIVERGENCE FIELDS (IN)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)                  ! SPECTRAL SCALAR FIELDS (IN)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)                 ! GRID VECTOR FIELDS     (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFVOR(:),YDFDIV(:)             ! GRID VECTOR FIELDS :VORTICITY AND DIVERGENCE     (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)                    ! GRID SCALAR FIELDS     (OUT)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU_EW(:),YDFV_EW(:)             ! GRID VECTOR FIELDS DERIVATIVES EW (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR_NS(:), YDFSCALAR_EW(:)  ! GRID SCALAR FIELDS DERIVATIVES EW AND NS (OUT)

INTEGER(KIND=JPIM),   INTENT(IN), OPTIONAL  :: KRESOL
PROCEDURE(FSPGL_INTF), POINTER, INTENT(IN), OPTIONAL  :: FSPGL_PROC

END SUBROUTINE INV_TRANS_FIELD_API

END INTERFACE
