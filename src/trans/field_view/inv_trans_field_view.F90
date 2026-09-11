! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE INV_TRANS_FIELD_VIEW(KRESOL, &
                               & YDSPSCALAR, YDSPVOR, YDSPDIV, &
                               & YDGPSCALAR, YDGPU, YDGPV,     &
                               & YDGPVOR,YDGPDIV, &
                               & YDGPSCALAR_NS, YDGPSCALAR_EW, YDGPU_EW, YDGPV_EW, &
                               & FSPGL_PROC)

!**** *INV_TRANS_FIELD_VIEW* - Field view interface to inverse spectral transform

!     Purpose.
!     --------
!        Allow to call INV_TRANS with a list of fields from field API

!**   Interface.
!     ----------
!     CALL INV_TRANS_FIELD_VIEW(...)

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

USE ISO_FORTRAN_ENV, ONLY : INT32
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW

USE YOMHOOK, ONLY : LHOOK,   DR_HOOK, JPHOOK
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW_GET_DATA_PTR, FIELD_VIEW_GET_IVSET_PTR
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW, LS_COUNT, LG_COUNT, LS, LG, IVSET_PTR, &
                                               & GET_NPROMA, GET_NFLD, GET_NSPEC2, GET_NBLK
USE TPM_DISTR, ONLY : DISTR_RESOL
USE PARKIND1, ONLY : JPRB, JPIM
USE ABORT_TRANS_MOD, ONLY : ABORT_TRANS
USE INV_TRANS_VIEW_CTL_MOD, ONLY: INV_TRANS_VIEW_CTL

IMPLICIT NONE

#include "fspgl_intf.h"
! INPUT
INTEGER(KIND=INT32),INTENT(IN) :: KRESOL
TYPE(FIELD_VIEW),INTENT(IN)  :: YDSPSCALAR(:)                  ! SPECTRAL SCALAR FIELDS (IN)
TYPE(FIELD_VIEW),INTENT(IN)  :: YDSPVOR(:), YDSPDIV(:)        ! SPECTRAL VECTOR FIELDS : VORTICITY AND DIVERGENCE FIELDS (IN)

! OUTPUT (note that views have intent(IN), i.e. the view itself is immutable, not the data)
TYPE(FIELD_VIEW),INTENT(IN)  :: YDGPSCALAR(:)                    ! GRID SCALAR FIELDS     (OUT)
TYPE(FIELD_VIEW),INTENT(IN)  :: YDGPU(:),YDGPV(:)                 ! GRID VECTOR FIELDS     (OUT)
TYPE(FIELD_VIEW),INTENT(IN)  :: YDGPVOR(:),YDGPDIV(:)             ! GRID VECTOR FIELDS :VORTICITY AND DIVERGENCE     (OUT)

TYPE(FIELD_VIEW),INTENT(IN)  :: YDGPSCALAR_EW(:), YDGPSCALAR_NS(:)  ! GRID SCALAR FIELDS DERIVATIVES EW AND NS (OUT)
TYPE(FIELD_VIEW),INTENT(IN)  :: YDGPU_EW(:),YDGPV_EW(:)             ! GRID VECTOR FIELDS DERIVATIVES EW (OUT)

PROCEDURE (FSPGL_INTF), POINTER, OPTIONAL, INTENT(IN)  :: FSPGL_PROC

! Local variables

! Lists of SPEC_VIEW and GRID_VIEW: intermediate representation of fields to facilitate copy to temporary arrays

TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVVOR(:), YLSPVDIV(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE :: YLGVU(:),YLGVV(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVVOR(:),YLGVDIV(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE :: YLGVU_EW(:),YLGVV_EW(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVSCALAR_NS(:), YLGVSCALAR_EW(:)

! b-set for inv-trans

TYPE(IVSET_PTR), ALLOCATABLE :: IVSETUV_LIST(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETSC_LIST(:)

INTEGER(KIND=JPIM)          :: IFLDSPVOR
INTEGER(KIND=JPIM)          :: JFLD, IFLD                             ! FIELD COUNTERS
LOGICAL                     :: LLSCDERS                               ! INDICATING IF DERIVATIVES OF SCALAR VARIABLES ARE REQ.
LOGICAL                     :: LLVORGP                                ! INDICATING IF GRID-POINT VORTICITY IS REQ.
LOGICAL                     :: LLDIVGP                                ! INDICATING IF GRID-POINT DIVERGENCE IS REQ.
LOGICAL                     :: LLUVDER                                ! INDICATING IF E-W DERIVATIVES OF U AND V ARE REQ.
INTEGER(KIND=JPIM)          :: NPROMA, NBLK
REAL(KIND=JPHOOK)           :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_VIEW',0,ZHOOK_HANDLE)

NPROMA              = GET_NPROMA(YDGPU, YDGPV, YDGPSCALAR)
NBLK                = GET_NBLK(YDGPU, YDGPV, YDGPSCALAR)

IFLDSPVOR= 0
JFLD= 0

LLSCDERS  = .FALSE.
LLVORGP = .FALSE.
LLDIVGP = .FALSE.
LLUVDER = .FALSE.

! 1. Vector fields transformation to grid space

! Preliminary checks
! Do we have vector fields?
IF (SIZE(YDGPU) > 0) THEN

  IF ((SIZE(YDGPU)/= SIZE(YDGPV)).OR.(SIZE(YDGPU)/= SIZE(YDSPDIV)).OR.(SIZE(YDGPU)/= SIZE(YDSPVOR))) THEN
    CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] The vector arrays have inconsistent sizes: YDGPU, YDGPV, YDSPDIV, YDSPVOR")
  ENDIF

  ! Convert list of spectral vector fields into a list of 2d FIELD_VIEW

  IFLDSPVOR = LS_COUNT(YDSPVOR)
  ALLOCATE(YLSPVVOR(IFLDSPVOR))
  ALLOCATE(YLSPVDIV(IFLDSPVOR))
  CALL LS(YDSPVOR, YLSPVVOR)
  CALL LS(YDSPDIV, YLSPVDIV)

  ! Convert list of grid-point vector fields into a list of 2d FIELD_VIEW
  ALLOCATE(YLGVU(LG_COUNT(YDGPU)))
  ALLOCATE(YLGVV(LG_COUNT(YDGPV)))
  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
    CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] inconsistent number of field_view for vectors")
  ENDIF

  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields.
  ALLOCATE(IVSETUV_LIST(SIZE(YDSPVOR)))
  DO JFLD=1,SIZE(YDSPVOR)
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPVOR(JFLD), IVSETUV_LIST(JFLD)%PTR)
  END DO
  ! Initialize b-set for vector fields data
  CALL LG(YDGPU, YLGVU, IVSETUV_LIST)
  CALL LG(YDGPV, YLGVV, IVSETUV_LIST)

  ! Output derivatives of vector fields
  IF (SIZE(YDGPU_EW) > 0 .AND. SIZE(YDGPV_EW) > 0)    THEN
    ALLOCATE(YLGVU_EW(LG_COUNT(YDGPU_EW)))
    ALLOCATE(YLGVV_EW(LG_COUNT(YDGPV_EW)))
    CALL LG(YDGPU_EW, YLGVU_EW, IVSETUV_LIST)
    CALL LG(YDGPV_EW, YLGVV_EW, IVSETUV_LIST)
 ENDIF

  ! Output divergence of vector fields
  IF (SIZE(YDGPDIV)  > 0) THEN
    ALLOCATE(YLGVDIV(LG_COUNT(YDGPDIV)))
    CALL LG(YDGPDIV, YLGVDIV, IVSETUV_LIST)
  ENDIF

  ! Output vorticity of vector fields
  IF (SIZE(YDGPVOR) > 0) THEN
    ALLOCATE(YLGVVOR(LG_COUNT(YDGPVOR)))
    CALL LG(YDGPVOR, YLGVVOR, IVSETUV_LIST)
  ENDIF
ENDIF

! 2. scalar fields transformation

IF (SIZE(YDSPSCALAR) > 0) THEN

  IF ((SIZE(YDSPSCALAR)/= SIZE(YDGPSCALAR))) CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] Inconsistent size &
                                                        & for YDSPSCALAR and YDGPSCALAR")

  ! Convert list of spectral scalar fields of any dimension into a list of 2d fields

! For LG we need the ivset of each grid point field,
! so we extract a matching list from the spectral fields
  ALLOCATE(IVSETSC_LIST(SIZE(YDGPSCALAR)))
  IFLD = 1
  DO JFLD=1,SIZE(YDSPSCALAR)
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPSCALAR(JFLD), IVSETSC_LIST(IFLD)%PTR)
    IVSETSC_LIST(IFLD)%BASE = YDSPSCALAR(JFLD)%IVSET_BASE
    IFLD = IFLD + 1
  END DO

  ALLOCATE(YLGVSCALAR(LG_COUNT(YDGPSCALAR)))
  CALL LG(YDGPSCALAR, YLGVSCALAR, IVSETSC_LIST)

  ALLOCATE(YLSPVSCALAR(LS_COUNT(YDSPSCALAR)))
  CALL LS(YDSPSCALAR, YLSPVSCALAR)

  IF (SIZE(YDGPSCALAR_NS) > 0 .AND. SIZE(YDGPSCALAR_EW) > 0) THEN
    ALLOCATE(YLGVSCALAR_NS(LG_COUNT(YDGPSCALAR_NS)))
    ALLOCATE(YLGVSCALAR_EW(LG_COUNT(YDGPSCALAR_EW)))
    CALL LG(YDGPSCALAR_NS, YLGVSCALAR_NS, IVSETSC_LIST)
    CALL LG(YDGPSCALAR_EW, YLGVSCALAR_EW, IVSETSC_LIST)
  ENDIF
ENDIF

! 3. CALL INV_TRANS_VIEW_CTL

IF (.NOT. ALLOCATED(YLSPVVOR))       ALLOCATE(YLSPVVOR(0))
IF (.NOT. ALLOCATED(YLSPVDIV))       ALLOCATE(YLSPVDIV(0))
IF (.NOT. ALLOCATED(YLSPVSCALAR))    ALLOCATE(YLSPVSCALAR(0))
IF (.NOT. ALLOCATED(YLGVU))          ALLOCATE(YLGVU(0))
IF (.NOT. ALLOCATED(YLGVV))          ALLOCATE(YLGVV(0))
IF (.NOT. ALLOCATED(YLGVVOR))        ALLOCATE(YLGVVOR(0))
IF (.NOT. ALLOCATED(YLGVDIV))        ALLOCATE(YLGVDIV(0))
IF (.NOT. ALLOCATED(YLGVSCALAR))     ALLOCATE(YLGVSCALAR(0))
IF (.NOT. ALLOCATED(YLGVU_EW))       ALLOCATE(YLGVU_EW(0))
IF (.NOT. ALLOCATED(YLGVV_EW))       ALLOCATE(YLGVV_EW(0))
IF (.NOT. ALLOCATED(YLGVSCALAR_NS))  ALLOCATE(YLGVSCALAR_NS(0))
IF (.NOT. ALLOCATED(YLGVSCALAR_EW))  ALLOCATE(YLGVSCALAR_EW(0))

CALL INV_TRANS_VIEW_CTL(NPROMA, NBLK, &
  & YLSPVVOR, YLSPVDIV, YLSPVSCALAR,&
  & YLGVU,YLGVV,&
  & YLGVVOR,YLGVDIV,&
  & YLGVSCALAR,&
  & YLGVU_EW,YLGVV_EW,&
  & YLGVSCALAR_NS, YLGVSCALAR_EW,&
  & FSPGL_PROC)

! 5. Final cleanup

IF (ALLOCATED(IVSETUV_LIST))  DEALLOCATE(IVSETUV_LIST)
IF (ALLOCATED(IVSETSC_LIST))  DEALLOCATE(IVSETSC_LIST)

! delete FIELD_VIEWS
IF (ALLOCATED(YLSPVVOR))    DEALLOCATE(YLSPVVOR)
IF (ALLOCATED(YLSPVDIV))    DEALLOCATE(YLSPVDIV)
IF (ALLOCATED(YLSPVSCALAR)) DEALLOCATE(YLSPVSCALAR)
IF (ALLOCATED(YLGVU))       DEALLOCATE(YLGVU)
IF (ALLOCATED(YLGVV))       DEALLOCATE(YLGVV)
IF (ALLOCATED(YLGVSCALAR))  DEALLOCATE(YLGVSCALAR)

IF (ALLOCATED(YLGVVOR))        DEALLOCATE(YLGVVOR)
IF (ALLOCATED(YLGVDIV))        DEALLOCATE(YLGVDIV)
IF (ALLOCATED(YLGVU_EW))       DEALLOCATE(YLGVU_EW)
IF (ALLOCATED(YLGVV_EW))       DEALLOCATE(YLGVV_EW)
IF (ALLOCATED(YLGVSCALAR_NS))  DEALLOCATE(YLGVSCALAR_NS)
IF (ALLOCATED(YLGVSCALAR_EW))  DEALLOCATE(YLGVSCALAR_EW)

IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_VIEW',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_FIELD_VIEW
