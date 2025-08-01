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


SUBROUTINE INV_TRANS_FIELD_API(YDFSPVOR,YDFSPDIV,YDFSPSCALAR, &
    & YDFU, YDFV, YDFVOR,YDFDIV,YDFSCALAR, &
    & YDFU_NS, YDFV_NS, YDFSCALAR_NS, YDFSCALAR_EW,& 
    & KSPEC, KPROMA, KGPBLKS, KGPTOT, KFLEVG, KFLEVL,KPROC,&
    & LDACC, &
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
!       YDFSPVOR(:)    - List of spectral vector fields (vorticity) 
!       YDFSPDIV(:)    - List of spectral vector fields (divergence)
!       YDFSPSCALAR(:) - List of spectral scalar fields 
!       KSPEC          - Number of spectral coefficients
!       KPROMA         - Blocking factor
!       KGPBLKS        - Number of blocks
!       KGPTOT         - Number of total grid points
!       KFLEVG         - Number of levels
!       KFLEVL         - Number of local levels
!       KPROC          - Processor ID
!       LDACC          - Field data on device
!       FSPGL_PROC     - procedure to be executed in fourier space
!                        before transposition

!      output
!       YDFU(:)        - List of grid-point vector fields (u)
!       YDFV(:)        - List of grid-point vector fields (v)
!       YDFVOR(:)      - List of grid-point vector fields (vorticity)
!       YDFDIV(:)      - List of grid-point vector fields (divergence)
!       YDFSCALAR(:)   - List of grid-point scalar fields
!       YDFU_NS(:)      - List of grid-point vector fields derivatives N-S (u)
!       YDFV_NS(:)      - List of grid-point vector fields derivatives N-S (v)
!       YDFSCALAR_NS(:) - List of grid-point scalar fields derivatives N-S
!       YDFSCALAR_EW(:) - List of grid-point scalar fields derivatives E-W

USE YOMHOOK, ONLY : LHOOK,   DR_HOOK, JPHOOK
USE FIELD_API_BASIC_TYPE_MOD, ONLY: FIELD_BASIC_PTR
USE PARKIND1, ONLY : JPIM, JPRB

#include "fspgl_intf.h"

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)        ! SPECTRAL VECTOR FIELDS : VORTICITY AND DIVERGENCE FIELDS (IN)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)                  ! SPECTRAL SCALAR FIELDS (IN)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)                 ! GRID VECTOR FIELDS     (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFVOR(:),YDFDIV(:)             ! GRID VECTOR FIELDS :VORTICITY AND DIVERGENCE     (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)                    ! GRID SCALAR FIELDS     (OUT)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU_NS(:),YDFV_NS(:)             ! GRID VECTOR FIELDS DERIVATIVES EW (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR_NS(:), YDFSCALAR_EW(:)  ! GRID SCALAR FIELDS DERIVATIVES EW AND NS (OUT)

INTEGER(KIND=JPIM),   INTENT(IN)            :: KSPEC
INTEGER(KIND=JPIM),   INTENT(IN)            :: KPROMA
INTEGER(KIND=JPIM),   INTENT(IN)            :: KGPBLKS
INTEGER(KIND=JPIM),   INTENT(IN)            :: KGPTOT
INTEGER(KIND=JPIM),   INTENT(IN)            :: KFLEVG
INTEGER(KIND=JPIM),   INTENT(IN)            :: KPROC
INTEGER(KIND=JPIM),   INTENT(IN)            :: KFLEVL
LOGICAL,              INTENT(IN), OPTIONAL  :: LDACC
PROCEDURE (FSPGL_INTF),           OPTIONAL  :: FSPGL_PROC

END SUBROUTINE INV_TRANS_FIELD_API

END INTERFACE
