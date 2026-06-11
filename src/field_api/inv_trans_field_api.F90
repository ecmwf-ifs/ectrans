! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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
USE ECTRANS_FIELD_API_MOD, ONLY : FIELD_GRID, FIELD_SPEC, SPEC_VIEW, GRID_VIEW, LS_COUNT, LG_COUNT, LS, LG, &
                                & GET_LAYOUT_S, GET_LAYOUT_G, IVSET_PTR
USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

#include "fspgl_intf.h"

INTEGER(KIND=JPIM),   INTENT(IN), OPTIONAL  :: KRESOL
TYPE(FIELD_SPEC),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)                  ! SPECTRAL SCALAR FIELDS (IN)
TYPE(FIELD_SPEC),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)        ! SPECTRAL VECTOR FIELDS : VORTICITY AND DIVERGENCE FIELDS (IN)

TYPE(FIELD_GRID),INTENT(INOUT), OPTIONAL  :: YDFSCALAR(:)                    ! GRID SCALAR FIELDS     (OUT)
TYPE(FIELD_GRID),INTENT(INOUT), OPTIONAL  :: YDFU(:),YDFV(:)                 ! GRID VECTOR FIELDS     (OUT)
TYPE(FIELD_GRID),INTENT(INOUT), OPTIONAL  :: YDFVOR(:),YDFDIV(:)             ! GRID VECTOR FIELDS :VORTICITY AND DIVERGENCE     (OUT)

TYPE(FIELD_GRID),INTENT(INOUT), OPTIONAL  :: YDFSCALAR_NS(:), YDFSCALAR_EW(:)  ! GRID SCALAR FIELDS DERIVATIVES EW AND NS (OUT)
TYPE(FIELD_GRID),INTENT(INOUT), OPTIONAL  :: YDFU_EW(:),YDFV_EW(:)             ! GRID VECTOR FIELDS DERIVATIVES EW (OUT)

INTEGER(KIND=JPIM),   INTENT(IN)            :: KGPTOT

PROCEDURE (FSPGL_INTF), POINTER, INTENT(IN), OPTIONAL  :: FSPGL_PROC

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
INTEGER(KIND=JPIM)          :: C
LOGICAL                     :: LLSCDERS                               ! INDICATING IF DERIVATIVES OF SCALAR VARIABLES ARE REQ.
LOGICAL                     :: LLVORGP                                ! INDICATING IF GRID-POINT VORTICITY IS REQ.
LOGICAL                     :: LLDIVGP                                ! INDICATING IF GRID-POINT DIVERGENCE IS REQ.
LOGICAL                     :: LLUVDER                                ! INDICATING IF E-W DERIVATIVES OF U AND V ARE REQ.
LOGICAL                     :: LDACC                                  ! INDICATING IF WE ARE RUNNING ON THE GPU
REAL(KIND=JPHOOK)           :: ZHOOK_HANDLE
INTEGER(KIND=JPIM)          :: NSPEC2
INTEGER(KIND=JPIM)          :: NBLK
INTEGER(KIND=JPIM)          :: NPROMA
INTEGER(KIND=JPIM)          :: NFLEVG

#include "inv_trans.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_API',0,ZHOOK_HANDLE)

NBLK = 0 
NPROMA = 0 
NSPEC2 = 0

IF (PRESENT(YDFU)) CALL GET_LAYOUT_G(YDFU, NBLK, NPROMA)
IF (PRESENT(YDFV)) CALL GET_LAYOUT_G(YDFV, NBLK, NPROMA)
IF (PRESENT(YDFSCALAR)) CALL GET_LAYOUT_G(YDFSCALAR,NBLK, NPROMA)
IF (PRESENT(YDFSPVOR)) CALL GET_LAYOUT_S(YDFSPVOR, NSPEC2)
IF (PRESENT(YDFSPDIV)) CALL GET_LAYOUT_S(YDFSPDIV, NSPEC2)
IF (PRESENT(YDFSPSCALAR)) CALL GET_LAYOUT_S(YDFSPSCALAR, NSPEC2)
IF (PRESENT(YDFSPVOR)) CALL GET_LAYOUT_S(YDFSPVOR, NSPEC2)

! Assert that NBLK, NPROMA and NSPEC2 are properly valued after getting the field(s) layout
IF (.NOT.(NBLK > 0))    CALL ABOR1("[INV_TRANS_FIELD_API] NBLK must be strictly positive. &
                                        &  One or more field arguments might be missing.")
IF (.NOT.(NPROMA > 0))  CALL ABOR1("[INV_TRANS_FIELD_API] NPROMA must be strictly positive. &
                                        &  One or more field arguments might be missing.")
IF (.NOT.(NSPEC2 > 0))  CALL ABOR1("[INV_TRANS_FIELD_API] NSPEC2 must be strictly positive. &
                                        &  One or more field arguments might be missing.")

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

! We are still relying on INV_TRANS, which require to have all the data are on CPU.
! So we force using data on the host 

LDACC = .FALSE.

! 1. Vector fields transformation to grid space

! Preliminary checks

IF (PRESENT(YDFU)      .NEQV. PRESENT(YDFV))     CALL ABOR1("[INV_TRANS_FIELD_API]  YDFU and YDFV &
                                                           & must be provided together")
IF (PRESENT(YDFSPDIV)  .NEQV. PRESENT(YDFSPVOR)) CALL ABOR1("[INV_TRANS_FIELD_API]  YDFSPDIV and YDFSPVOR &
                                                           & must be provided together")
IF (PRESENT(YDFU) .AND. .NOT. PRESENT(YDFSPVOR)) CALL ABOR1("[INV_TRANS_FIELD_API] YDFU and YDFSPVOR &
                                                           & must be provided together")
IF (PRESENT(YDFU) .AND. .NOT. PRESENT(YDFSPDIV)) CALL ABOR1("[INV_TRANS_FIELD_API] YDFU and YDFSPDIV &
                                                           & must be provided together")

! Do we have vector fields?
IF (PRESENT(YDFU)) THEN

  IF ((SIZE(YDFU)/= SIZE(YDFV)).OR.(SIZE(YDFU)/= SIZE(YDFSPDIV)).OR.(SIZE(YDFU)/= SIZE(YDFSPVOR))) THEN
    CALL ABOR1("[INV_TRANS_FIELD_API] The vector arrays have inconsistent sizes: YDFU, YDFV, YDFSPDIV, YDFSPVOR")
  ENDIF

  ! Convert list of spectral vector fields into a list of 2d FIELD_VIEW

  IFLDSPVOR = LS_COUNT(YDFSPVOR)
  ALLOCATE(YLSPVVOR(IFLDSPVOR))
  ALLOCATE(YLSPVDIV(IFLDSPVOR))
  ! Convert list of grid-point vector fields into a list of 2d FIELD_VIEW
  ALLOCATE(YLGVU(LG_COUNT(YDFU)))
  ALLOCATE(YLGVV(LG_COUNT(YDFV)))
  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
    CALL ABOR1("[INV_TRANS_FIELD_API] inconsistent number of field_view for vectors")
  ENDIF

  NFLEVG = SIZE (YLGVU) / SIZE (YDFU)
  IUVG = SIZE(YDFU)

  LLUVDER  = .FALSE.
  LLVORGP = .FALSE.
  LLDIVGP = .FALSE.
  LLSCDERS = .FALSE.

  IUVDIM = 2

  ! Output derivatives of vector fields
  IF (PRESENT(YDFU_EW) .AND. PRESENT(YDFV_EW))    THEN
    LLUVDER = .TRUE.
    IUVDIM = IUVDIM + 2
    ALLOCATE(YLGVU_EW(LG_COUNT(YDFU_EW)))
    ALLOCATE(YLGVV_EW(LG_COUNT(YDFV_EW)))
 ENDIF

  ! Output divergence of vector fields
  IF (PRESENT(YDFDIV)) THEN
    LLDIVGP = .TRUE.
    IUVDIM = IUVDIM + 1
    ALLOCATE(YLGVDIV(LG_COUNT(YDFDIV)))
  ENDIF

  ! Output vorticity of vector fields
  IF (PRESENT(YDFVOR)) THEN
    LLVORGP = .TRUE.
    IUVDIM = IUVDIM + 1
    ALLOCATE(YLGVVOR(LG_COUNT(YDFVOR)))
  ENDIF

  ! allocate temporary vector field arrays in spectral space
  ALLOCATE(ZPSPVOR(IFLDSPVOR,NSPEC2))
  ALLOCATE(ZPSPDIV(IFLDSPVOR,NSPEC2))

  ! allocate temporary vector field array in grid space
  ALLOCATE(ZPGPUV(NPROMA,NFLEVG, IUVG * IUVDIM,NBLK))

  ! allocate 'b-set' for vector fields
  ALLOCATE(IVSETUV_LIST(SIZE(YDFSPVOR)))
  DO JFLD=1,SIZE(YDFSPVOR)
    IVSETUV_LIST(JFLD)%PTR => YDFSPVOR(JFLD)%IVSET
  END DO

  C = LS(YDFSPVOR, YLSPVVOR, LDACC, .TRUE.)
  C = LS(YDFSPDIV, YLSPVDIV, LDACC, .TRUE.)

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
  C = LG(YDFU, YLGVU, IVSETUV_LIST, LDACC, .TRUE.)
  ALLOCATE(IVSETUV(NFLEVG))
  DO JFLD=1,IUVG
    DO JLEV=1,NFLEVG
      ID = JLEV + (JFLD -1) * NFLEVG
      IF (JFLD .EQ. 1) THEN
        IVSETUV(JLEV) = YLGVU(ID)%IVSET
      ENDIF
      IF (IVSETUV(JLEV) .NE. YLGVU(ID)%IVSET) CALL ABOR1("[INV_TRANS_FIELD_API] ivsetuv inconsistent with ylgvu%ivset")
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

IF (PRESENT(YDFSPSCALAR) .NEQV. PRESENT(YDFSCALAR)) CALL ABOR1("[INV_TRANS_FIELD_API]  YDFSPSCALAR and YDFSCALAR &
                                                              & must be provided together")

IF (PRESENT(YDFSPSCALAR)) THEN

  IF ((SIZE(YDFSPSCALAR)/= SIZE(YDFSCALAR))) CALL ABOR1("[INV_TRANS_FIELD_API] Inconsistent size &
                                                        & for YDFSPSCALAR and YDFSCALAR")

  ! Convert list of spectral scalar fields of any domension into a list of 2d fields
  IFLDSPSC = LS_COUNT(YDFSPSCALAR)
  ALLOCATE(YLSPVSCALAR(IFLDSPSC))

  ALLOCATE(YLGVSCALAR(LG_COUNT(YDFSCALAR)))

  IFLDXG = SIZE(YLGVSCALAR) ! NUMBER OF OUTPUT SCALAR FIELDS IN GRID SPACE
  ! count the number of fields present on the processor
  C = LS(YDFSPSCALAR, YLSPVSCALAR, LDACC,.TRUE.)
  IFLDXL = 0
  DO JFLD = 1, IFLDSPSC
    IF (ASSOCIATED(YLSPVSCALAR(JFLD)%P)) THEN
      IFLDXL = IFLDXL + 1
    ENDIF
  END DO
  ISCDIM = 1
  IF (PRESENT(YDFSCALAR_NS) .AND. PRESENT(YDFSCALAR_EW)) THEN
    LLSCDERS = .TRUE.
    ISCDIM = ISCDIM + 2
    ALLOCATE(YLGVSCALAR_NS(LG_COUNT(YDFSCALAR_NS)))
    ALLOCATE(YLGVSCALAR_EW(LG_COUNT(YDFSCALAR_EW)))
 ENDIF

! Allocate scalar field array in spectral space
  ALLOCATE(ZPSPSC2(IFLDXL,NSPEC2))

! Allocate scalar field array in grid space
  ALLOCATE(ZPGP2(NPROMA,IFLDXG * ISCDIM,NBLK))

! allocate 'b-set' for scalar fields
  ALLOCATE(IVSETSC_LIST(SIZE(YDFSCALAR)))
  IFLD = 1
  DO JFLD=1,SIZE(YDFSPSCALAR)
      IVSETSC_LIST(IFLD)%PTR => YDFSPSCALAR(JFLD)%IVSET
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
  C = LG(YDFSCALAR, YLGVSCALAR, IVSETSC_LIST, LDACC,.TRUE.)
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

! 4. Copy back temporary array data into grid-point fields

! remove garbage at the end of arrays
IEND = KGPTOT - NPROMA * (NBLK - 1)

IF (IUVG>0) ZPGPUV (IEND+1:, :, :, NBLK) = 0
IF (IFLDXG>0)  ZPGP2 (IEND+1:, :, NBLK) = 0

! copy vector fields

IF (IUVG>0) THEN

  IOFFSET = 0
  ! copy vorticity
  IF (LLVORGP) THEN
      C = LG(YDFVOR, YLGVVOR, IVSETUV_LIST, LDACC, .FALSE.)
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
      C = LG(YDFDIV, YLGVDIV, IVSETUV_LIST, LDACC, .FALSE.)
      DO JFLD=1,IUVG
        DO JLEV=1,NFLEVG
          ID = JLEV + (JFLD -1) * NFLEVG
          ZZ2_1=>YLGVDIV(ID)%P
          ZZ2_1(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
          ENDDO
      ENDDO

    IOFFSET = IOFFSET + 1
  ENDIF

  ! copy u and v
  C = LG(YDFU, YLGVU, IVSETUV_LIST, LDACC, .FALSE.)
  C = LG(YDFV, YLGVV, IVSETUV_LIST, LDACC, .FALSE.)


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
    C = LG(YDFU_EW, YLGVU_EW, IVSETUV_LIST, LDACC, .FALSE.)
    C = LG(YDFV_EW, YLGVV_EW, IVSETUV_LIST, LDACC, .FALSE.)

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
    C = LG(YDFSCALAR, YLGVSCALAR, IVSETSC_LIST, LDACC,.FALSE.)
    DO JFLD=1, IFLDXG
      ZZ2_1=>YLGVSCALAR(JFLD)%P(:,:)
      ZZ2_1(:,:) = ZPGP2(:,JFLD,:)
    ENDDO

  ! copy spectral scalar fields derivatives

  IF (LLSCDERS) THEN
    C = LG(YDFSCALAR_NS, YLGVSCALAR_NS, IVSETSC_LIST, LDACC, .FALSE.)
    C = LG(YDFSCALAR_EW, YLGVSCALAR_EW, IVSETSC_LIST, LDACC, .FALSE.)

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

IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_API',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_FIELD_API
