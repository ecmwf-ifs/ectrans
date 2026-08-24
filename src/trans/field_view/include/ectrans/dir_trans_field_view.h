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

SUBROUTINE DIR_TRANS_FIELD_VIEW(KRESOL,                       &
                              & YDGPSCALAR, YDGPU, YDGPV,     &
                              & YDSPSCALAR, YDSPVOR, YDSPDIV)

!**** *DIR_TRANS_FIELD_VIEW* - Field view interface to direct spectral transform
!
!     Purpose.
!     --------
!        Allow to call DIR_TRANS with a list of FIELD_VIEW fields
!
!     Explicit arguments :
!     --------------------
!      input
!       KRESOL         - The resolution identifier
!       YDGPSCALAR(:)  - List of grid-point scalar fields
!       YDGPU(:)       - List of grid-point vector fields (u)
!       YDGPV(:)       - List of grid-point vector fields (v)
!
!      output
!       YDSPSCALAR(:)  - List of spectral scalar fields
!       YDSPVOR(:)     - List of spectral vector fields (vorticity)
!       YDSPDIV(:)     - List of spectral vector fields (divergence)

USE ISO_FORTRAN_ENV, ONLY : INT32
USE ECTRANS_FIELD_VIEW_MOD, ONLY : FIELD_VIEW

IMPLICIT NONE

INTEGER(KIND=INT32), INTENT(IN)    :: KRESOL
TYPE(FIELD_VIEW),   INTENT(IN)    :: YDGPSCALAR(:), YDGPU(:), YDGPV(:)
TYPE(FIELD_VIEW),   INTENT(INOUT) :: YDSPSCALAR(:), YDSPVOR(:), YDSPDIV(:)

END SUBROUTINE DIR_TRANS_FIELD_VIEW

END INTERFACE
