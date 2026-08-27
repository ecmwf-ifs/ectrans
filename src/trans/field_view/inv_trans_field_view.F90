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

IMPLICIT NONE

#include "fspgl_intf.h"
INTEGER(KIND=INT32),   INTENT(IN) :: KRESOL
TYPE(FIELD_VIEW),INTENT(IN)  :: YDSPSCALAR(:)                  ! SPECTRAL SCALAR FIELDS (IN)
TYPE(FIELD_VIEW),INTENT(IN)  :: YDSPVOR(:), YDSPDIV(:)        ! SPECTRAL VECTOR FIELDS : VORTICITY AND DIVERGENCE FIELDS (IN)

TYPE(FIELD_VIEW),INTENT(INOUT)  :: YDGPSCALAR(:)                    ! GRID SCALAR FIELDS     (OUT)
TYPE(FIELD_VIEW),INTENT(INOUT)  :: YDGPU(:),YDGPV(:)                 ! GRID VECTOR FIELDS     (OUT)
TYPE(FIELD_VIEW),INTENT(INOUT)  :: YDGPVOR(:),YDGPDIV(:)             ! GRID VECTOR FIELDS :VORTICITY AND DIVERGENCE     (OUT)

TYPE(FIELD_VIEW),INTENT(INOUT)  :: YDGPSCALAR_EW(:), YDGPSCALAR_NS(:)  ! GRID SCALAR FIELDS DERIVATIVES EW AND NS (OUT)
TYPE(FIELD_VIEW),INTENT(INOUT)  :: YDGPU_EW(:),YDGPV_EW(:)             ! GRID VECTOR FIELDS DERIVATIVES EW (OUT)

PROCEDURE (FSPGL_INTF), POINTER, OPTIONAL, INTENT(IN)  :: FSPGL_PROC

! Local variables

LOGICAL :: LLFSPGL_PROC

! List of FIELD_VIEW: intermediate representation of fields to facilitate copy to temporary arrays

TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVVOR(:), YLSPVDIV(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE :: YLGVU(:),YLGVV(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVVOR(:),YLGVDIV(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE :: YLGVU_EW(:),YLGVV_EW(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVSCALAR_NS(:), YLGVSCALAR_EW(:)

! Temporary arrays for inv_trans
REAL(KIND=JPRB), POINTER :: ZPSPVOR(:,:),ZPSPDIV(:,:)  ! spectral vector fields (in)
REAL(KIND=JPRB), POINTER :: ZPSPSC2(:,:)               ! spectral scalar fields (in)
REAL(KIND=JPRB), POINTER :: ZPGPUV(:,:,:,:)            ! grid vector fields (out)
REAL(KIND=JPRB), POINTER :: ZPGP2(:,:,:)               ! grid scalar fields (out)

REAL(KIND=JPRB), POINTER :: ZZ1_1(:)
REAL(KIND=JPRB), POINTER :: ZZ1_2(:)
REAL(KIND=JPRB), POINTER :: ZZ2_1(:,:)
REAL(KIND=JPRB), POINTER :: ZZ2_2(:,:)

! b-set for inv-trans
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETUV(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETSC2(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETUV_LIST(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETSC_LIST(:)

INTEGER(KIND=JPIM)          :: IFLDXG
INTEGER(KIND=JPIM)          :: IFLDXL
INTEGER(KIND=JPIM)          :: IFLDSPVOR
INTEGER(KIND=JPIM)          :: IFLDSPSC
INTEGER(KIND=JPIM)          :: IUVG
INTEGER(KIND=JPIM)          :: ISCDIM
INTEGER(KIND=JPIM)          :: IUVDIM
INTEGER(KIND=JPIM)          :: ID,IOFFSET,JLEV
INTEGER(KIND=JPIM)          :: IEND
INTEGER(KIND=JPIM)          :: JFLD, IFLD                             ! FIELD COUNTERS
LOGICAL                     :: LLSCDERS                               ! INDICATING IF DERIVATIVES OF SCALAR VARIABLES ARE REQ.
LOGICAL                     :: LLVORGP                                ! INDICATING IF GRID-POINT VORTICITY IS REQ.
LOGICAL                     :: LLDIVGP                                ! INDICATING IF GRID-POINT DIVERGENCE IS REQ.
LOGICAL                     :: LLUVDER                                ! INDICATING IF E-W DERIVATIVES OF U AND V ARE REQ.
INTEGER(KIND=JPIM)          :: NFLEVG
INTEGER(KIND=JPIM)          :: NGPTOT

INTEGER(KIND=JPIM) :: NFIELD_UV, NFIELD_SCALAR
INTEGER(KIND=JPIM) :: NPROMA, NBLK, NFIELD_TOTAL_UV, NSPEC2

REAL(KIND=JPHOOK)           :: ZHOOK_HANDLE

#include "inv_trans.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_VIEW',0,ZHOOK_HANDLE)


NFIELD_UV           = SIZE(YDGPU)
NFIELD_SCALAR       = SIZE(YDGPSCALAR)
NPROMA              = GET_NPROMA(YDGPU, YDGPV, YDGPSCALAR)
NBLK                = GET_NBLK(YDGPU, YDGPV, YDGPSCALAR)
NFIELD_TOTAL_UV     = GET_NFLD(YDGPU)
NSPEC2              = GET_NSPEC2(YDSPVOR, YDSPDIV, YDSPSCALAR)

IFLDXG= 0
IFLDXL= 0
IFLDSPVOR= 0
IFLDSPSC= 0
IUVG  = 0
ISCDIM = 0
IUVDIM = 0
ID= 0
IOFFSET= 0
JLEV= 0
JFLD= 0
IEND= 0
NFLEVG = 0
LLSCDERS  = .FALSE.
LLVORGP = .FALSE.
LLDIVGP = .FALSE.
LLUVDER = .FALSE.

LLFSPGL_PROC = .FALSE.
IF (PRESENT(FSPGL_PROC)) THEN
  IF (ASSOCIATED(FSPGL_PROC)) THEN
    LLFSPGL_PROC = .TRUE.
  ENDIF
ENDIF

IF ((SIZE(YDGPU) /= SIZE(YDGPV)).OR.(SIZE(YDGPU) /= SIZE(YDSPDIV)).OR.(SIZE(YDGPU) /= SIZE(YDSPVOR))) THEN
  CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] The vector arrays have inconsistent sizes: YDGPU, YDGPV, YDSPDIV, YDSPVOR")
ENDIF
IF (SIZE(YDGPVOR) > 0 .AND. SIZE(YDGPVOR) /= SIZE(YDSPVOR)) THEN
  CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] YDGPVOR and YDSPVOR must have the same size")
ENDIF
IF (SIZE(YDGPDIV) > 0 .AND. SIZE(YDGPDIV) /= SIZE(YDSPDIV)) THEN
  CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] YDGPDIV and YDSPDIV must have the same size")
ENDIF
IF (SIZE(YDGPU_EW) > 0 .OR. SIZE(YDGPV_EW) > 0) THEN
  IF (SIZE(YDGPU_EW) /= SIZE(YDGPV_EW) .OR. SIZE(YDGPU_EW) /= SIZE(YDGPU)) THEN
    CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] YDGPU_EW, YDGPV_EW and YDGPU must have the same size")
  ENDIF
ENDIF
IF ((SIZE(YDSPSCALAR)/= SIZE(YDGPSCALAR))) THEN
  CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] Inconsistent size for YDSPSCALAR and YDGPSCALAR")
ENDIF



! 1. Vector fields transformation to grid space

! Do we have vector fields?
IF (SIZE(YDGPU) > 0) THEN

  ! Convert list of spectral vector fields into a list of 2d FIELD_VIEW

  IFLDSPVOR = LS_COUNT(YDSPVOR)
  ALLOCATE(YLSPVVOR(IFLDSPVOR))
  ALLOCATE(YLSPVDIV(IFLDSPVOR))
  ! Convert list of grid-point vector fields into a list of 2d FIELD_VIEW
  ALLOCATE(YLGVU(LG_COUNT(YDGPU)))
  ALLOCATE(YLGVV(LG_COUNT(YDGPV)))
  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
    CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] inconsistent number of field_view for vectors")
  ENDIF

  NFLEVG = SIZE (YLGVU) / SIZE (YDGPU)
  IUVG = SIZE(YDGPU)

  LLUVDER  = .FALSE.
  LLVORGP = .FALSE.
  LLDIVGP = .FALSE.
  LLSCDERS = .FALSE.

  IUVDIM = 2

  ! Output derivatives of vector fields
  IF (SIZE(YDGPU_EW) > 0 .AND. SIZE(YDGPV_EW) > 0)    THEN
    LLUVDER = .TRUE.
    IUVDIM = IUVDIM + 2
    ALLOCATE(YLGVU_EW(LG_COUNT(YDGPU_EW)))
    ALLOCATE(YLGVV_EW(LG_COUNT(YDGPV_EW)))
  ENDIF

  ! Output divergence of vector fields
  IF (SIZE(YDGPDIV)  > 0) THEN
    LLDIVGP = .TRUE.
    IUVDIM = IUVDIM + 1
    ALLOCATE(YLGVDIV(LG_COUNT(YDGPDIV)))
  ENDIF

  ! Output vorticity of vector fields
  IF (SIZE(YDGPVOR) > 0) THEN
    LLVORGP = .TRUE.
    IUVDIM = IUVDIM + 1
    IF (.NOT. LLDIVGP) IUVDIM = IUVDIM + 1
    ALLOCATE(YLGVVOR(LG_COUNT(YDGPVOR)))
  ENDIF

  ! allocate temporary vector field arrays in spectral space
  ALLOCATE(ZPSPVOR(IFLDSPVOR,NSPEC2))
  ALLOCATE(ZPSPDIV(IFLDSPVOR,NSPEC2))

  ! allocate temporary vector field array in grid space
  ALLOCATE(ZPGPUV(NPROMA,NFLEVG, IUVG * IUVDIM,NBLK))

  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields.
  ALLOCATE(IVSETUV_LIST(SIZE(YDSPVOR)))
  DO JFLD=1,SIZE(YDSPVOR)
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPVOR(JFLD), IVSETUV_LIST(JFLD)%PTR)
  END DO

  CALL LS(YDSPVOR, YLSPVVOR)
  CALL LS(YDSPDIV, YLSPVDIV)

  ! Copy list of 2d views of spectral vector fields into temporary arrays
  DO JFLD=1,IFLDSPVOR
    IF (ASSOCIATED(YLSPVVOR(JFLD)%P)) THEN
      ZZ1_1=>YLSPVVOR(JFLD)%P
      ZZ1_2=>YLSPVDIV(JFLD)%P
      ZPSPVOR(JFLD,:) = ZZ1_1(:)
      ZPSPDIV(JFLD,:) = ZZ1_2(:)
    ENDIF
  ENDDO

  ! Initialize b-set for vector fields data
  CALL LG(YDGPU, YLGVU, IVSETUV_LIST)
  ALLOCATE(IVSETUV(NFLEVG))
  DO JFLD=1,IUVG
    DO JLEV=1,NFLEVG
      ID = JLEV + (JFLD - 1) * NFLEVG
      IF (JFLD .EQ. 1) THEN
        IVSETUV(JLEV) = YLGVU(ID)%IVSET
      ENDIF
      IF (IVSETUV(JLEV) .NE. YLGVU(ID)%IVSET) THEN
        CALL ABORT_TRANS("[INV_TRANS_FIELD_VIEW] ivsetuv inconsistent with ylgvu%ivset")
      ENDIF
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

IF (SIZE(YDSPSCALAR) > 0) THEN

  ! Convert list of spectral scalar fields of any domension into a list of 2d fields
  IFLDSPSC = LS_COUNT(YDSPSCALAR)
  ALLOCATE(YLSPVSCALAR(IFLDSPSC))

  ALLOCATE(YLGVSCALAR(LG_COUNT(YDGPSCALAR)))

  IFLDXG = SIZE(YLGVSCALAR) ! NUMBER OF OUTPUT SCALAR FIELDS IN GRID SPACE
  ! count the number of fields present on the processor
  CALL LS(YDSPSCALAR, YLSPVSCALAR)
  IFLDXL = 0
  DO JFLD = 1, IFLDSPSC
    IF (ASSOCIATED(YLSPVSCALAR(JFLD)%P)) THEN
      IFLDXL = IFLDXL + 1
    ENDIF
  END DO
  ISCDIM = 1
  IF (SIZE(YDGPSCALAR_NS) > 0 .AND. SIZE(YDGPSCALAR_EW) > 0) THEN
    LLSCDERS = .TRUE.
    ISCDIM = ISCDIM + 2
    ALLOCATE(YLGVSCALAR_NS(LG_COUNT(YDGPSCALAR_NS)))
    ALLOCATE(YLGVSCALAR_EW(LG_COUNT(YDGPSCALAR_EW)))
  ENDIF

  ! Allocate scalar field array in spectral space
  ALLOCATE(ZPSPSC2(IFLDXL,NSPEC2))

  ! Allocate scalar field array in grid space
  ALLOCATE(ZPGP2(NPROMA,IFLDXG * ISCDIM,NBLK))

  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields
  ALLOCATE(IVSETSC_LIST(SIZE(YDGPSCALAR)))
  IFLD = 1
  DO JFLD=1,SIZE(YDSPSCALAR)
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPSCALAR(JFLD), IVSETSC_LIST(IFLD)%PTR)
    IFLD = IFLD + 1
  END DO

  ! Copy list of of spectral scalar fields into temporary arrays (1d copy thanks to field_view)
  ID = 1
  DO JFLD = 1,IFLDSPSC
    IF (ASSOCIATED(YLSPVSCALAR(JFLD)%P)) THEN
      ZZ1_1=>YLSPVSCALAR(JFLD)%P
      ZPSPSC2(ID,:) = ZZ1_1(:)
      ID = ID + 1
    ENDIF
  ENDDO

  ! compute ´b-set´ for scalar-fields
  CALL LG(YDGPSCALAR, YLGVSCALAR, IVSETSC_LIST)
  ALLOCATE(IVSETSC2(IFLDXG))
  DO JFLD=1, IFLDXG
    IVSETSC2(JFLD) = YLGVSCALAR(JFLD)%IVSET
  ENDDO
ELSE
  !No scalar field provided
  IFLDXG = 0
  ZPGP2=>NULL()
  ZPSPSC2=>NULL()
ENDIF

! 3. CALL INV_TRANS  using the regular interface and the temporary arrays

! We have to perform separated calls for nvfortran
IF (ASSOCIATED(ZPGP2) .AND. ASSOCIATED(ZPGPUV)) THEN
  IF (LLFSPGL_PROC) THEN
    CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                  & PSPSC2=ZPSPSC2,PGP2=ZPGP2,KVSETSC2=IVSETSC2, &
                  & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                  & KPROMA=NPROMA, FSPGL_PROC=FSPGL_PROC, KRESOL=KRESOL)
  ELSE
    CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                  & PSPSC2=ZPSPSC2,PGP2=ZPGP2, KVSETSC2=IVSETSC2, &
                  & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                  & KPROMA=NPROMA, KRESOL=KRESOL)
  ENDIF
ELSE IF (ASSOCIATED(ZPGP2)) THEN
  IF (LLFSPGL_PROC) THEN
    CALL INV_TRANS(PSPSC2=ZPSPSC2,PGP2=ZPGP2,KVSETSC2=IVSETSC2, &
                  & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                  & KPROMA=NPROMA, FSPGL_PROC=FSPGL_PROC, KRESOL=KRESOL)
  ELSE
    CALL INV_TRANS(PSPSC2=ZPSPSC2,PGP2=ZPGP2, KVSETSC2=IVSETSC2, &
                  & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                  & KPROMA=NPROMA, KRESOL=KRESOL)
  ENDIF
ELSE IF (ASSOCIATED(ZPGPUV)) THEN
  IF (LLFSPGL_PROC) THEN
    CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                  & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                  & KPROMA=NPROMA, FSPGL_PROC=FSPGL_PROC, KRESOL=KRESOL)
  ELSE
    CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                  & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                  & KPROMA=NPROMA, KRESOL=KRESOL)
  ENDIF
ENDIF

! Get NGPTOT from TPM_DISTR module's DIST_RESOL(KRESOL)
NGPTOT = DISTR_RESOL(KRESOL)%NGPTOT

! 4. Copy back temporary array data into grid-point fields

! remove garbage at the end of arrays
IEND = NGPTOT - NPROMA * (NBLK - 1)

IF (IUVG>0) ZPGPUV (IEND+1:, :, :, NBLK) = 0
IF (IFLDXG>0)  ZPGP2 (IEND+1:, :, NBLK) = 0

! copy vector fields

IF (IUVG>0) THEN

  IOFFSET = 0
  ! copy vorticity
  IF (LLVORGP) THEN
    CALL LG(YDGPVOR, YLGVVOR, IVSETUV_LIST)
    DO JFLD=1,IUVG
      DO JLEV=1,NFLEVG
        ID = JLEV + (JFLD -1) * NFLEVG
        ZZ2_1=>YLGVVOR(ID)%P
        ZZ2_1(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
      ENDDO
    ENDDO

    IOFFSET = IOFFSET + 1
  ENDIF

  ! copy divergence
  IF (LLDIVGP) THEN
    CALL LG(YDGPDIV, YLGVDIV, IVSETUV_LIST)
    DO JFLD=1,IUVG
      DO JLEV=1,NFLEVG
        ID = JLEV + (JFLD -1) * NFLEVG
        ZZ2_1=>YLGVDIV(ID)%P
        ZZ2_1(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
      ENDDO
    ENDDO
  ENDIF
  IF (LLDIVGP .OR. LLVORGP) IOFFSET = IOFFSET + 1

  ! copy u and v
  CALL LG(YDGPU, YLGVU, IVSETUV_LIST)
  CALL LG(YDGPV, YLGVV, IVSETUV_LIST)


  DO JFLD=1,IUVG
    DO JLEV=1,NFLEVG
      ID = JLEV + (JFLD -1) * NFLEVG
      ZZ2_1=>YLGVU(ID)%P
      ZZ2_2=>YLGVV(ID)%P
      ZZ2_1(:,:) =  ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:)
      ZZ2_2(:,:) =  ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:)
    ENDDO
  ENDDO

  IOFFSET = IOFFSET + 2

  ! copy u and v derivatives
  IF (LLUVDER) THEN
    CALL LG(YDGPU_EW, YLGVU_EW, IVSETUV_LIST)
    CALL LG(YDGPV_EW, YLGVV_EW, IVSETUV_LIST)

    DO JFLD=1,IUVG
      DO JLEV=1,NFLEVG
        ID = JLEV + (JFLD -1) * NFLEVG
        ZZ2_1=>YLGVU_EW(ID)%P
        ZZ2_2=>YLGVV_EW(ID)%P
        ZZ2_1(:,:) =  ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:)
        ZZ2_2(:,:) =  ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:)
      ENDDO
    ENDDO
  ENDIF
ENDIF

IF (IFLDXG > 0) THEN
  ! copy spectral scalar fields
  CALL LG(YDGPSCALAR, YLGVSCALAR, IVSETSC_LIST)
  DO JFLD=1, IFLDXG
    ZZ2_1=>YLGVSCALAR(JFLD)%P(:,:)
    ZZ2_1(:,:) = ZPGP2(:,JFLD,:)
  ENDDO

  ! copy spectral scalar fields derivatives
  IF (LLSCDERS) THEN
    CALL LG(YDGPSCALAR_NS, YLGVSCALAR_NS, IVSETSC_LIST)
    CALL LG(YDGPSCALAR_EW, YLGVSCALAR_EW, IVSETSC_LIST)

    DO JFLD=1,IFLDXG
      ZZ2_1=>YLGVSCALAR_NS(JFLD)%P
      ZZ2_2=>YLGVSCALAR_EW(JFLD)%P
      ZZ2_1(:,:) = ZPGP2(:, JFLD+IFLDXG,:)
      ZZ2_2(:,:) = ZPGP2(:, JFLD+(2*IFLDXG),:)
    ENDDO
  ENDIF
ENDIF

! 5. Final cleanup

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

IF (ALLOCATED(YLGVVOR))        DEALLOCATE(YLGVVOR)
IF (ALLOCATED(YLGVDIV))        DEALLOCATE(YLGVDIV)
IF (ALLOCATED(YLGVU_EW))       DEALLOCATE(YLGVU_EW)
IF (ALLOCATED(YLGVV_EW))       DEALLOCATE(YLGVV_EW)
IF (ALLOCATED(YLGVSCALAR_NS))  DEALLOCATE(YLGVSCALAR_NS)
IF (ALLOCATED(YLGVSCALAR_EW))  DEALLOCATE(YLGVSCALAR_EW)

IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_VIEW',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_FIELD_VIEW
