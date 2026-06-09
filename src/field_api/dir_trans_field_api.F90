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
                               & GET_LAYOUT_S, GET_LAYOUT_G, IVSET_PTR
USE PARKIND1  ,ONLY : JPIM,JPRB, JPRD

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KRESOL
TYPE(FIELD_GRID),INTENT(IN), OPTIONAL  :: YDFSCALAR(:), YDFU(:), YDFV(:)
TYPE(FIELD_SPEC),INTENT(INOUT), OPTIONAL  :: YDFSPSCALAR(:), YDFSPVOR(:), YDFSPDIV(:)

! Local variables

! List of FIELD_VIEW: intermediate representation of fields to facilitate copy to temporary arrays
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVVOR(:), YLSPVDIV(:)
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVU(:),YLGVV(:)
TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVSCALAR(:)

! Temporary arrays for dir_trans
REAL(KIND=JPRB),POINTER :: ZPSPVOR(:,:),ZPSPDIV(:,:)  ! spectral vector fields (out)
REAL(KIND=JPRB),POINTER :: ZPSPSC2(:,:)               ! spectral scalar fields(out)
REAL(KIND=JPRB),POINTER :: ZPGPUV(:,:,:,:)            ! grid vector fields (in)
REAL(KIND=JPRB),POINTER :: ZPGP2(:,:,:)               ! grid scalar fields (in)

REAL(KIND=JPRB), POINTER :: ZZ1_1(:)
REAL(KIND=JPRB), POINTER :: ZZ1_2(:)
REAL(KIND=JPRB), POINTER :: ZZ2_1(:,:)
REAL(KIND=JPRB), POINTER :: ZZ2_2(:,:)
REAL(KIND=JPRB):: S

! b-set for dir-trans
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETUV(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETSC2(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETUV_LIST(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETSC_LIST(:)

INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NPROMA
INTEGER(KIND=JPIM) :: NBLK
INTEGER(KIND=JPIM) :: KGPTOT
INTEGER(KIND=JPIM) :: NFLEVG

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
INTEGER(KIND=JPIM) :: C
LOGICAL            :: LDACC     ! INDICATING IF WE ARE RUNNING ON THE GPU
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

#include "dir_trans.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',0,ZHOOK_HANDLE)

NBLK = 0 
NPROMA = 0 
NSPEC2 = 0

IF (PRESENT(YDFU))        CALL GET_LAYOUT_G(YDFU, NBLK, NPROMA)
IF (PRESENT(YDFV))        CALL GET_LAYOUT_G(YDFV, NBLK, NPROMA)
IF (PRESENT(YDFSCALAR))   CALL GET_LAYOUT_G(YDFSCALAR,NBLK, NPROMA)
IF (PRESENT(YDFSPVOR))    CALL GET_LAYOUT_S(YDFSPVOR, NSPEC2)
IF (PRESENT(YDFSPDIV))    CALL GET_LAYOUT_S(YDFSPDIV, NSPEC2)
IF (PRESENT(YDFSPSCALAR)) CALL GET_LAYOUT_S(YDFSPSCALAR, NSPEC2)
IF (PRESENT(YDFSPVOR))    CALL GET_LAYOUT_S(YDFSPVOR, NSPEC2)

! Assert that NBLK, NPROMA and NSPEC2 are properly valued after getting the field(s) layout
IF (.NOT.(NBLK > 0))    CALL ABOR1("[DIR_TRANS_FIELD_API] NBLK must be strictly positive. &
                                        &  One or more field arguments might be missing.")
IF (.NOT.(NPROMA > 0))  CALL ABOR1("[DIR_TRANS_FIELD_API] NPROMA must be strictly positive. &
                                        &  One or more field arguments might be missing.")
IF (.NOT.(NSPEC2 > 0))  CALL ABOR1("[DIR_TRANS_FIELD_API] NSPEC2 must be strictly positive. &
                                        &  One or more field arguments might be missing.")


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

! We are still relying on DIR_TRANS, which require to have all the data are on CPU.
! So we force using data on the host
LDACC = .FALSE.

! 1. Vector fields transformation to spectral space

! Do we have vector fields?
IF (PRESENT(YDFU)) THEN

  IF (.NOT. PRESENT(YDFV))     CALL ABOR1("[DIR_TRANS_FIELD_API] YDFU and YDFV must be provided together for vector transform")
  IF (.NOT. PRESENT(YDFSPVOR)) CALL ABOR1("[DIR_TRANS_FIELD_API] YDFSPVOR must be provided for vector transform")
  IF (.NOT. PRESENT(YDFSPDIV)) CALL ABOR1("[DIR_TRANS_FIELD_API] YDFSPDIV must be provided for vector transform")

  IF ((SIZE(YDFU)/= SIZE(YDFV)).OR.(SIZE(YDFU)/= SIZE(YDFSPDIV)).OR.(SIZE(YDFU)/= SIZE(YDFSPVOR))) THEN
     CALL ABOR1("[DIR_TRANS_FIELD_API] The vector arrays have inconsistent sizes: YDFU, YDFV, YDFSPDIV, YDFSPVOR")
  ENDIF

  ! Convert list of spectral vector fields into a list of 2d FIELD_VIEW
  IFLDSPVOR = LS_COUNT(YDFSPVOR)

  ALLOCATE(YLSPVVOR(IFLDSPVOR))
  ALLOCATE(YLSPVDIV(IFLDSPVOR))

  ! Convert list of grid-point vector fields into a list of 2d FIELD_VIEW
  ALLOCATE(YLGVU(LG_COUNT(YDFU)))
  ALLOCATE(YLGVV(LG_COUNT(YDFV)))
  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
     CALL ABOR1("[DIR_TRANS_FIELD_API] inconsistent number of field_view for vectors")
  ENDIF
  NFLEVG = SIZE (YLGVU) / SIZE (YDFU)
  IUVG = SIZE(YDFU)

  IUVDIM = 2

  ! allocate temporary vector field arrays in spectral space
  ALLOCATE(ZPSPVOR(IFLDSPVOR,NSPEC2))
  ALLOCATE(ZPSPDIV(IFLDSPVOR,NSPEC2))

  ! allocate temporary vector field array in grid space
  ALLOCATE(ZPGPUV(NPROMA,NFLEVG, IUVG * IUVDIM,NBLK))

  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields.
  ALLOCATE(IVSETUV_LIST(IUVG))
  DO JFLD=1,IUVG
    IVSETUV_LIST(JFLD)%PTR => YDFSPVOR(JFLD)%IVSET
  ENDDO

  C = LG(YDFU, YLGVU, IVSETUV_LIST, LDACC, .TRUE.)
  C = LG(YDFV, YLGVV, IVSETUV_LIST, LDACC, .TRUE.)

  ! Copy list of 2d views of grid point vector fields into temporary arrays
  IOFFSET = 0
  DO JFLD=1,IUVG
    DO JLEV=1,NFLEVG
      ID = JLEV + (JFLD -1) * NFLEVG
      ZZ2_1=>YLGVU(ID)%P
      ZZ2_2=>YLGVV(ID)%P

      ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:) = ZZ2_1(:,:)
      ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:) = ZZ2_2(:,:)
    ENDDO
  ENDDO

  ALLOCATE(IVSETUV(NFLEVG))
  DO JFLD=1,IUVG
    DO JLEV=1,NFLEVG
      ID = JLEV + (JFLD -1) * NFLEVG
      IF (JFLD .EQ. 1) THEN
        IVSETUV(JLEV) = YLGVU(ID)%IVSET
      ENDIF
      IF (IVSETUV(JLEV) .NE. YLGVV(ID)%IVSET)  CALL ABOR1("[DIR_TRANS_FIELD_API] ivsetuv inconsistent with ylgvv%ivset")
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
IF (PRESENT(YDFSPSCALAR) .NEQV. PRESENT(YDFSCALAR))  CALL ABOR1("[DIR_TRANS_FIELD_API] YDFSPSCALAR and YDFSCALAR &
                                                               & must be provided together")

! Do we have scalar fields?
IF (PRESENT(YDFSPSCALAR)) THEN
  IF ((SIZE(YDFSPSCALAR)/= SIZE(YDFSCALAR)))  CALL ABOR1("[DIR_TRANS_FIELD_API] Inconsistent size &
                                                         & for YDFSPSCALAR and YDFSCALAR")

  ! Convert list of spectral scalar fields of any dimension into a list of 2d fields
  ALLOCATE(YLGVSCALAR(LG_COUNT(YDFSCALAR)))

  IFLDXG = SIZE(YLGVSCALAR)

  IFLDSPSC = LS_COUNT(YDFSPSCALAR)
  ALLOCATE(YLSPVSCALAR(IFLDSPSC))

  ! count the number of fields present on the processor
  C = LS(YDFSPSCALAR, YLSPVSCALAR, LDACC,.TRUE.)
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
  ALLOCATE(IVSETSC_LIST(SIZE(YDFSPSCALAR)))
  DO JFLD=1,SIZE(YDFSPSCALAR)
    IVSETSC_LIST(JFLD)%PTR => YDFSPSCALAR(JFLD)%IVSET
  ENDDO

  ! Copy list of scalar fields into temporary arrays (2d copy thanks to field_view)
  C = LG(YDFSCALAR, YLGVSCALAR, IVSETSC_LIST, LDACC,.TRUE.)
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

    C = LS(YDFSPVOR, YLSPVVOR, LDACC, .FALSE.)
    C = LS(YDFSPDIV, YLSPVDIV, LDACC, .FALSE.)

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

   C = LS(YDFSPSCALAR, YLSPVSCALAR, LDACC,.FALSE.)
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

! delete FIELD_VIEWS
IF (ALLOCATED(YLSPVVOR))    DEALLOCATE(YLSPVVOR)
IF (ALLOCATED(YLSPVDIV))    DEALLOCATE(YLSPVDIV)
IF (ALLOCATED(YLSPVSCALAR)) DEALLOCATE(YLSPVSCALAR)
IF (ALLOCATED(YLGVU))       DEALLOCATE(YLGVU)
IF (ALLOCATED(YLGVV))       DEALLOCATE(YLGVV)
IF (ALLOCATED(YLGVSCALAR))  DEALLOCATE(YLGVSCALAR)

IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_FIELD_API
