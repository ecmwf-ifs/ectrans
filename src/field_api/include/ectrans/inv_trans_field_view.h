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
                               & YDGPVOR, YDGPDIV,             &
                               & YDGPSCALAR_NS, YDGPSCALAR_EW, YDGPU_EW, YDGPV_EW, &
                               & FSPGL_PROC)

!**** *INV_TRANS_FIELD_VIEW* - Field API interface to inverse spectral transform
!
!     Purpose.
!     --------
!        Allow to call INV_TRANS with a list of FIELD_VIEW fields
!
!     Explicit arguments :
!     --------------------
!      input
!       KRESOL           The resolution identifier
!       YDSPSCALAR(:)  - List of spectral scalar fields
!       YDSPVOR(:)     - List of spectral vector fields (vorticity)
!       YDSPDIV(:)     - List of spectral vector fields (divergence)
!       FSPGL_PROC     - procedure to be executed in fourier space
!                        before transposition
!
!      output
!       YDGPSCALAR(:)    - List of grid-point scalar fields
!       YDGPU(:)         - List of grid-point vector fields (u)
!       YDGPV(:)         - List of grid-point vector fields (v)
!       YDGPVOR(:)       - List of grid-point vector fields (vorticity)
!       YDGPDIV(:)       - List of grid-point vector fields (divergence)
!       YDGPSCALAR_NS(:) - List of grid-point scalar fields derivatives N-S
!       YDGPSCALAR_EW(:) - List of grid-point scalar fields derivatives E-W
!       YDGPU_EW(:)      - List of grid-point vector fields derivatives E-W (u)
!       YDGPV_EW(:)      - List of grid-point vector fields derivatives E-W (v)

USE PARKIND1, ONLY : JPIM
USE ECTRANS_FIELD_VIEW_MOD, ONLY : FIELD_VIEW

IMPLICIT NONE

#include "fspgl_intf.h"

INTEGER(KIND=JPIM), INTENT(IN)    :: KRESOL
TYPE(FIELD_VIEW),   INTENT(IN)    :: YDSPSCALAR(:)
TYPE(FIELD_VIEW),   INTENT(IN)    :: YDSPVOR(:), YDSPDIV(:)

TYPE(FIELD_VIEW),   INTENT(INOUT) :: YDGPSCALAR(:)
TYPE(FIELD_VIEW),   INTENT(INOUT) :: YDGPU(:), YDGPV(:)
TYPE(FIELD_VIEW),   INTENT(INOUT) :: YDGPVOR(:), YDGPDIV(:)

TYPE(FIELD_VIEW),   INTENT(INOUT) :: YDGPSCALAR_NS(:), YDGPSCALAR_EW(:)
TYPE(FIELD_VIEW),   INTENT(INOUT) :: YDGPU_EW(:), YDGPV_EW(:)

PROCEDURE(FSPGL_INTF), POINTER, OPTIONAL, INTENT(IN) :: FSPGL_PROC

END SUBROUTINE INV_TRANS_FIELD_VIEW

END INTERFACE
