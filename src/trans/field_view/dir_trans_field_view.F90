! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DIR_TRANS_FIELD_VIEW(KRESOL,                       &
                              & YDGPSCALAR, YDGPU, YDGPV,     &
                              & YDSPSCALAR, YDSPVOR, YDSPDIV)

USE ISO_FORTRAN_ENV, ONLY : INT32
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW


USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARKIND1 ,ONLY : JPRB, JPIM
USE ABORT_TRANS_MOD, ONLY : ABORT_TRANS
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW_GET_IVSET_PTR
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY: SPEC_VIEW, GRID_VIEW, LS_COUNT, LG_COUNT, LS, LG, IVSET_PTR, &
                                              & GET_NPROMA, GET_NFLD, GET_NSPEC2, GET_NBLK

IMPLICIT NONE

! Arguments
INTEGER(KIND=INT32), INTENT(IN)      :: KRESOL
TYPE(FIELD_VIEW),   INTENT(IN)      :: YDGPSCALAR(:), YDGPU(:), YDGPV(:)
TYPE(FIELD_VIEW),   INTENT(INOUT)   :: YDSPSCALAR(:), YDSPVOR(:), YDSPDIV(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: NFIELD_UV, NFIELD_SCALAR
INTEGER(KIND=JPIM) :: NPROMA, NBLK, NFIELD_TOTAL_UV, NSPEC2

! Local variables

! List of FIELD_VIEW: intermediate representation of fields to facilitate copy to temporary arrays
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVSCALAR(:)
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVVOR(:), YLSPVDIV(:)

TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVSCALAR(:)
TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVU(:),YLGVV(:)

! Temporary arrays for dir_trans
REAL(KIND=JPRB),POINTER :: ZPSPVOR(:,:),ZPSPDIV(:,:)  ! spectral vector fields (out)
REAL(KIND=JPRB),POINTER :: ZPSPSC2(:,:)               ! spectral scalar fields(out)
REAL(KIND=JPRB),POINTER :: ZPGPUV(:,:,:,:)            ! grid vector fields (in)
REAL(KIND=JPRB),POINTER :: ZPGP2(:,:,:)               ! grid scalar fields (in)

REAL(KIND=JPRB), POINTER :: ZZ1_1(:)
REAL(KIND=JPRB), POINTER :: ZZ1_2(:)
REAL(KIND=JPRB), POINTER :: ZZ2_1(:,:)
REAL(KIND=JPRB), POINTER :: ZZ2_2(:,:)

! b-set for dir-trans
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETUV(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETSC2(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETUV_LIST(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETSC_LIST(:)

INTEGER(KIND=JPIM) :: IFLDXG
INTEGER(KIND=JPIM) :: IFLDXL
INTEGER(KIND=JPIM) :: IFLDSPVOR
INTEGER(KIND=JPIM) :: IFLDSPSC
INTEGER(KIND=JPIM) :: IUVG
INTEGER(KIND=JPIM) :: IUVDIM
INTEGER(KIND=JPIM) :: ID
INTEGER(KIND=JPIM) :: IOFFSET
INTEGER(KIND=JPIM) :: JLEV      ! Level counter
INTEGER(KIND=JPIM) :: JFLD      ! Field counter

INTEGER(KIND=JPIM) :: KFLEVG
#include "dir_trans.h"

IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_VIEW',0,ZHOOK_HANDLE)

NFIELD_UV           = SIZE(YDGPU)
NFIELD_SCALAR       = SIZE(YDGPSCALAR) 
NPROMA              = GET_NPROMA(YDGPU, YDGPV, YDGPSCALAR)
NBLK                = GET_NBLK(YDGPU, YDGPV, YDGPSCALAR)
NFIELD_TOTAL_UV     = GET_NFLD(YDGPU)
NSPEC2              = GET_NSPEC2(YDSPVOR, YDSPDIV, YDSPSCALAR)

IFLDXG  = 0
IFLDXL = 0
IFLDSPVOR= 0
IFLDSPSC= 0
IUVG = 0
JFLD  = 0
IUVDIM = 0
ID = 0
IOFFSET = 0
JLEV = 0
JFLD = 0

! 1. Vector fields transformation to spectral space

! Preliminary checks

! Do we have vector fields?
IF (SIZE(YDGPU) > 0) THEN

  IF ((SIZE(YDGPU) /= SIZE(YDGPV)).OR.(SIZE(YDGPU) /= SIZE(YDSPDIV)).OR.(SIZE(YDGPU)/= SIZE(YDSPVOR))) THEN
    CALL ABORT_TRANS("[DIR_TRANS_FIELD_VIEW] The vector arrays have inconsistent sizes: YDGPU, YDGPV, YDSPDIV, YDSPVOR")
  ENDIF

  ! Convert list of spectral vector fields into a list of 2d FIELD_VIEW
  IFLDSPVOR = LS_COUNT(YDSPVOR)

  ALLOCATE(YLSPVVOR(IFLDSPVOR))
  ALLOCATE(YLSPVDIV(IFLDSPVOR))

  ! Convert list of grid-point vector fields into a list of 2d FIELD_VIEW
  ALLOCATE(YLGVU(LG_COUNT(YDGPU)))
  ALLOCATE(YLGVV(LG_COUNT(YDGPV)))
  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
    CALL ABORT_TRANS("[DIR_TRANS_FIELD_VIEW] inconsistent number of field_view for vectors")
  ENDIF
  KFLEVG = SIZE (YLGVU) / SIZE (YDGPU)
  IUVG = SIZE(YDGPU)

  IUVDIM = 2

  ! allocate temporary vector field arrays in spectral space
  ALLOCATE(ZPSPVOR(IFLDSPVOR,NSPEC2))
  ALLOCATE(ZPSPDIV(IFLDSPVOR,NSPEC2))

  ! allocate temporary vector field array in grid space
  ALLOCATE(ZPGPUV(NPROMA,KFLEVG, IUVG * IUVDIM,NBLK))

  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields.
  ALLOCATE(IVSETUV_LIST(IUVG))
  DO JFLD=1,IUVG
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPVOR(JFLD), IVSETUV_LIST(JFLD)%PTR)
  ENDDO

  CALL LG(YDGPU, YLGVU, IVSETUV_LIST)
  CALL LG(YDGPV, YLGVV, IVSETUV_LIST)

  ! Copy list of 2d views of grid point vector fields into temporary arrays
  IOFFSET = 0
  DO JFLD=1,IUVG
    DO JLEV=1,KFLEVG
      ID = JLEV + (JFLD -1) * KFLEVG
      ZZ2_1=>YLGVU(ID)%P
      ZZ2_2=>YLGVV(ID)%P

      ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:) = ZZ2_1(:,:)
      ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:) = ZZ2_2(:,:)
    ENDDO
  ENDDO

  ALLOCATE(IVSETUV(KFLEVG))
  DO JFLD=1,IUVG
    DO JLEV=1,KFLEVG
      ID = JLEV + (JFLD -1) * KFLEVG
      IF (JFLD .EQ. 1) THEN
        IVSETUV(JLEV) = YLGVU(ID)%IVSET
      ENDIF
      IF (IVSETUV(JLEV) .NE. YLGVV(ID)%IVSET)  CALL ABORT_TRANS("[DIR_TRANS_FIELD_VIEW] ivsetuv inconsistent with ylgvv%ivset")
    ENDDO
  ENDDO
ELSE
  ! No vector field provided
  IUVG = 0  
  ZPGPUV=>NULL()
  ZPSPVOR=>NULL()
  ZPSPDIV=>NULL()
ENDIF

! 2. scalar fields transformation

! Preliminary checks

! Do we have scalar fields?
IF (SIZE(YDSPSCALAR) > 0 ) THEN
  IF ((SIZE(YDSPSCALAR)/= SIZE(YDGPSCALAR)))  CALL ABORT_TRANS("[DIR_TRANS_FIELD_VIEW] Inconsistent size &
                                                               & for YDSPSCALAR and YDGPSCALAR")

  ! Convert list of spectral scalar fields of any dimension into a list of 2d fields
  ALLOCATE(YLGVSCALAR(LG_COUNT(YDGPSCALAR)))

  IFLDXG = SIZE(YLGVSCALAR)

  IFLDSPSC = LS_COUNT(YDSPSCALAR)
  ALLOCATE(YLSPVSCALAR(IFLDSPSC))

  ! count the number of fields present on the processor
  CALL LS(YDSPSCALAR, YLSPVSCALAR)
  IFLDXL = 0
  DO JFLD = 1, IFLDSPSC
    IF (ASSOCIATED(YLSPVSCALAR(JFLD)%P)) IFLDXL = IFLDXL + 1
  END DO
   ! Allocate temporary scalar field array in spectral space
  ALLOCATE(ZPSPSC2(IFLDXL,NSPEC2))

  ! Allocate temporary scalar field array in grid space
  ALLOCATE(ZPGP2(NPROMA,IFLDXG,NBLK))

  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields
  ALLOCATE(IVSETSC_LIST(SIZE(YDSPSCALAR)))
  DO JFLD=1,SIZE(YDSPSCALAR)
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPSCALAR(JFLD), IVSETSC_LIST(JFLD)%PTR)
  ENDDO

  ! Copy list of scalar fields into temporary arrays (2d copy thanks to field_view)
  CALL LG(YDGPSCALAR, YLGVSCALAR, IVSETSC_LIST)
  ALLOCATE(IVSETSC2(IFLDXG))
  DO JFLD=1, IFLDXG
    ZZ2_1=>YLGVSCALAR(JFLD)%P
    ZPGP2(:,JFLD,:) = ZZ2_1(:,:)
    IVSETSC2(JFLD) = YLGVSCALAR(JFLD)%IVSET
  ENDDO

ELSE
  !No scalar field provided  
  IFLDXG = 0
  ZPGP2=>NULL()
  ZPSPSC2=>NULL()
ENDIF

! 3. CALL DIR_TRANS using the regular interface and the temporary arrays

! We have to perform separated calls for nvfortran
IF (ASSOCIATED(ZPGP2) .AND. ASSOCIATED(ZPGPUV)) THEN
  CALL DIR_TRANS(PSPVOR = ZPSPVOR,PSPDIV = ZPSPDIV,PGPUV = ZPGPUV,KVSETUV = IVSETUV, &
                & PSPSC2 = ZPSPSC2,PGP2 = ZPGP2, KVSETSC2 = IVSETSC2, &
                & KPROMA = NPROMA, KRESOL = KRESOL)
ELSE IF (ASSOCIATED(ZPGP2)) THEN
  CALL DIR_TRANS(PSPSC2 = ZPSPSC2,PGP2 = ZPGP2, KVSETSC2 = IVSETSC2, &
                & KPROMA = NPROMA, KRESOL = KRESOL)
ELSE IF (ASSOCIATED(ZPGPUV)) THEN
  CALL DIR_TRANS(PSPVOR = ZPSPVOR,PSPDIV = ZPSPDIV,PGPUV = ZPGPUV,KVSETUV = IVSETUV, &
                & KPROMA = NPROMA, KRESOL = KRESOL)
ENDIF
! 4. Copy back temporary array data into spectral fields

! copy spectral vorticity and divergence
IF (IUVG>0) THEN
  CALL LS(YDSPVOR, YLSPVVOR)
  CALL LS(YDSPDIV, YLSPVDIV)

  DO JFLD=1,IFLDSPVOR
    IF (ASSOCIATED(YLSPVVOR(JFLD)%P)) THEN
      ZZ1_1=>YLSPVVOR(JFLD)%P
      ZZ1_2=>YLSPVDIV(JFLD)%P
      ZZ1_1(:) = ZPSPVOR(JFLD,:)
      ZZ1_2(:) = ZPSPDIV(JFLD,:)
    ENDIF
  ENDDO
ENDIF

! copy spectral scalar fields
IF (IFLDSPSC > 0) THEN

  CALL LS(YDSPSCALAR, YLSPVSCALAR)
  ID = 1
  DO JFLD = 1, IFLDSPSC
    IF (ASSOCIATED(YLSPVSCALAR(JFLD)%P)) THEN
      ZZ1_1=>YLSPVSCALAR(JFLD)%P
      ZZ1_1(:) = ZPSPSC2(ID,:)
      ID = ID + 1
    ENDIF
  ENDDO
ENDIF

! 5. Final cleanup

! delete temporary arrays

IF (ASSOCIATED(ZPSPVOR)) DEALLOCATE(ZPSPVOR)
IF (ASSOCIATED(ZPSPDIV)) DEALLOCATE(ZPSPDIV)
IF (ASSOCIATED(ZPSPSC2)) DEALLOCATE(ZPSPSC2)
IF (ASSOCIATED(ZPGPUV))  DEALLOCATE(ZPGPUV)
IF (ASSOCIATED(ZPGP2))   DEALLOCATE(ZPGP2)
IF (ALLOCATED(IVSETUV))  DEALLOCATE(IVSETUV)
IF (ALLOCATED(IVSETSC2)) DEALLOCATE(IVSETSC2)
IF (ALLOCATED(IVSETUV_LIST))  DEALLOCATE(IVSETUV_LIST)
IF (ALLOCATED(IVSETSC_LIST))  DEALLOCATE(IVSETSC_LIST)

! delete FIELD_VIEWS
IF (ALLOCATED(YLSPVVOR))    DEALLOCATE(YLSPVVOR)
IF (ALLOCATED(YLSPVDIV))    DEALLOCATE(YLSPVDIV)
IF (ALLOCATED(YLSPVSCALAR)) DEALLOCATE(YLSPVSCALAR)
IF (ALLOCATED(YLGVU))       DEALLOCATE(YLGVU)
IF (ALLOCATED(YLGVV))       DEALLOCATE(YLGVV)
IF (ALLOCATED(YLGVSCALAR))  DEALLOCATE(YLGVSCALAR)

IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_VIEW',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_FIELD_VIEW
