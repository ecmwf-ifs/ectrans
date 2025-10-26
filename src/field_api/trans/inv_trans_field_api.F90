! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE INV_TRANS_FIELD_API(YDFSPVOR,YDFSPDIV,YDFSPSCALAR, &
                             & YDFU, YDFV, YDFVOR,YDFDIV,YDFSCALAR, &
                             & YDFU_EW, YDFV_EW, YDFSCALAR_NS, YDFSCALAR_EW,&
                             & KSPEC, KPROMA, KGPBLKS, KGPTOT, KFLEVG, KFLEVL, KPROC,&
                             & LDACC,&
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
!       LDACC          - Field and temporary data on the device
!       FSPGL_PROC     - procedure to be executed in fourier space
!                        before transposition

!      output
!       YDFU(:)        - List of grid-point vector fields (u)
!       YDFV(:)        - List of grid-point vector fields (v)
!       YDFVOR(:)      - List of grid-point vector fields (vorticity)
!       YDFDIV(:)      - List of grid-point vector fields (divergence)
!       YDFSCALAR(:)   - List of grid-point scalar fields
!       YDFU_EW(:)      - List of grid-point vector fields derivatives E-W (u)
!       YDFV_EW(:)      - List of grid-point vector fields derivatives E-W (v)
!       YDFSCALAR_NS(:) - List of grid-point scalar fields derivatives N-S
!       YDFSCALAR_EW(:) - List of grid-point scalar fields derivatives E-W

USE YOMHOOK, ONLY : LHOOK,   DR_HOOK, JPHOOK
USE FIELD_API_BASIC_TYPE_MOD, ONLY: FIELD_BASIC_PTR
USE FIELD_API_ECTRANS_MOD
USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

#include "fspgl_intf.h"

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)        ! SPECTRAL VECTOR FIELDS : VORTICITY AND DIVERGENCE FIELDS (IN)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)                  ! SPECTRAL SCALAR FIELDS (IN)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)                 ! GRID VECTOR FIELDS     (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFVOR(:),YDFDIV(:)             ! GRID VECTOR FIELDS :VORTICITY AND DIVERGENCE     (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)                    ! GRID SCALAR FIELDS     (OUT)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU_EW(:),YDFV_EW(:)             ! GRID VECTOR FIELDS DERIVATIVES EW (OUT)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR_NS(:), YDFSCALAR_EW(:)  ! GRID SCALAR FIELDS DERIVATIVES EW AND NS (OUT)

INTEGER(KIND=JPIM),   INTENT(IN)            :: KSPEC
INTEGER(KIND=JPIM),   INTENT(IN)            :: KPROMA
INTEGER(KIND=JPIM),   INTENT(IN)            :: KGPBLKS
INTEGER(KIND=JPIM),   INTENT(IN)            :: KGPTOT
INTEGER(KIND=JPIM),   INTENT(IN)            :: KFLEVG
INTEGER(KIND=JPIM),   INTENT(IN)            :: KFLEVL
INTEGER(KIND=JPIM),   INTENT(IN)            :: KPROC
LOGICAL,              INTENT(IN), OPTIONAL  :: LDACC
PROCEDURE (FSPGL_INTF),           OPTIONAL  :: FSPGL_PROC

! Local variables

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

INTEGER(KIND=JPIM)          :: ISPUV
INTEGER(KIND=JPIM)          :: IFLDXG
INTEGER(KIND=JPIM)          :: IFLDXL
INTEGER(KIND=JPIM)          :: IFLDXGUV
INTEGER(KIND=JPIM)          :: IFLDXLUV
INTEGER(KIND=JPIM)          :: IFLDSPVOR
INTEGER(KIND=JPIM)          :: IFLDSPSC
INTEGER(KIND=JPIM)          :: IUVG
INTEGER(KIND=JPIM)          :: ISCDIM
INTEGER(KIND=JPIM)          :: IUVDIM
INTEGER(KIND=JPIM)          :: ID,IOFFSET,JLEV
INTEGER(KIND=JPIM)          :: IEND
INTEGER(KIND=JPIM)          :: JFLD                                   ! FIELD COUNTER
INTEGER(KIND=JPIM)          :: C
LOGICAL                     :: LLSCDERS                               ! INDICATING IF DERIVATIVES OF SCALAR VARIABLES ARE REQ.
LOGICAL                     :: LLVORGP                                ! INDICATING IF GRID-POINT VORTICITY IS REQ.
LOGICAL                     :: LLDIVGP                                ! INDICATING IF GRID-POINT DIVERGENCE IS REQ.
LOGICAL                     :: LLUVDER                                ! INDICATING IF E-W DERIVATIVES OF U AND V ARE REQ.
REAL(KIND=JPHOOK)           :: ZHOOK_HANDLE

#include "inv_trans.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_API',0,ZHOOK_HANDLE)

ISPUV = 0
IFLDXG= 0
IFLDXL= 0
IFLDXGUV= 0
IFLDXLUV= 0
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

! 1. Vector fields transformation to grid space

! Preliminary checks

IF (PRESENT(YDFU) .NEQV. PRESENT(YDFV)) CALL ABOR1("[INV_TRANS_FIELD_API]  YDFU and YDFV must be provided together")
IF (PRESENT(YDFSPDIV) .NEQV. PRESENT(YDFSPVOR)) CALL ABOR1("[INV_TRANS_FIELD_API]  YDFSPDIV and YDFSPVOR must be provided together")
IF (PRESENT(YDFU) .AND. .NOT. PRESENT(YDFSPVOR)) CALL ABOR1("[INV_TRANS_FIELD_API] YDFU and YDFSPVOR must be provided together")
IF (PRESENT(YDFU) .AND. .NOT. PRESENT(YDFSPDIV)) CALL ABOR1("[INV_TRANS_FIELD_API] YDFU and YDFSPDIV must be provided together")

! Do we have vector fields?
IF (PRESENT(YDFU)) THEN

  IF ((SIZE(YDFU)/= SIZE(YDFV)).OR.(SIZE(YDFU)/= SIZE(YDFSPDIV)).OR.(SIZE(YDFU)/= SIZE(YDFSPVOR))) THEN
    CALL ABOR1("[INV_TRANS_FIELD_API] The vector arrays have inconsitent sizes: YDFU, YDFV, YDFSPDIV, YDFSPVOR")
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
  IF (((SIZE (YLGVU) / SIZE (YDFU)) /= KFLEVG) .OR. ((SIZE (YLSPVVOR) / SIZE (YDFSPVOR)) /= KFLEVL)) THEN
    CALL ABOR1("[INV_TRANS_FIELD_API] inconsistent kflevg or kflevl")
  ENDIF

  IUVG = SIZE(YDFU)
  ISPUV = SIZE(YDFSPVOR)

  LLUVDER  = .FALSE.
  LLVORGP = .FALSE.
  LLDIVGP = .FALSE.
  LLSCDERS = .FALSE.

  IUVDIM = 2

  ! Output derivatives of vector fields
  IF (PRESENT(YDFU_EW) .AND. PRESENT(YDFV_EW))    THEN
    LLUVDER = .TRUE.
    IUVDIM = 5
    ALLOCATE(YLGVU_EW(LG_COUNT(YDFU_EW)))
    ALLOCATE(YLGVV_EW(LG_COUNT(YDFV_EW)))
 ENDIF

  ! Output divergence of vector fields
  IF (PRESENT(YDFDIV)) THEN
    LLDIVGP = .TRUE.
    IUVDIM = 5
    ALLOCATE(YLGVDIV(LG_COUNT(YDFDIV)))
  ENDIF

  ! Output vorticity of vector fields
  IF (PRESENT(YDFVOR)) THEN
    LLVORGP = .TRUE.
    IUVDIM = 6
    ALLOCATE(YLGVVOR(LG_COUNT(YDFVOR)))
  ENDIF

  ! allocate temporary vector field arrays in spectral space
  ALLOCATE(ZPSPVOR(IFLDSPVOR,KSPEC))
  ALLOCATE(ZPSPDIV(IFLDSPVOR,KSPEC))

  ! allocate temporary vector field array in grid space
  ALLOCATE(ZPGPUV(KPROMA,KFLEVG, IUVG * IUVDIM,KGPBLKS))

  ! allocate 'b-set' for vector fields
  ALLOCATE(IVSETUV(KFLEVG))

    ! temporary copies on gpu
  if ( LDACC ) THEN
    !$ACC ENTER DATA CREATE(ZPSPVOR,ZPSPDIV,ZPGPUV)
  ENDIF

    C = LS(YLSPVVOR, YDFSPVOR, LDACC, .TRUE.)
    C = LS(YLSPVDIV, YDFSPDIV, LDACC, .TRUE.)

    ! Copy list of 2d views of spectral vector fields into temporary arrays
    DO JFLD=1,IFLDSPVOR
    IF (ASSOCIATED(YLSPVVOR(JFLD)%P)) THEN
        ZZ1_1=>YLSPVVOR(JFLD)%P
        ZZ1_2=>YLSPVDIV(JFLD)%P
        IF (LDACC) THEN
          !$acc kernels present(ZPSPVOR,ZPSPDIV,ZZ1_1,ZZ1_2)
          ZPSPVOR(JFLD,:) = ZZ1_1(:)
          ZPSPDIV(JFLD,:) = ZZ1_2(:)
          !$acc end kernels
        ELSE
     	     ZPSPVOR(JFLD,:) = ZZ1_1(:)
           ZPSPDIV(JFLD,:) = ZZ1_2(:)
        ENDIF
      ENDIF
    ENDDO
  if (LDACC)  THEN
    !$ACC UPDATE SELF(ZPSPVOR,ZPSPDIV)
 endif
 
  ! Initialize b-set for vector fields data
  C = LG(YLGVU, YDFU, LDACC, .TRUE.)
  DO JFLD=1,IUVG
    DO JLEV=1,KFLEVG
     ID = JLEV + (JFLD -1) * KFLEVG
     IF (JFLD .EQ. 1) IVSETUV(JLEV) = YLGVU(ID)%IVSET
     IF (IVSETUV(JLEV) .NE. YLGVU(ID)%IVSET) CALL ABOR1("[INV_TRANS_FIELD_API] ivsetuv inconsistent with ylgvu%ivset")
    ENDDO
  ENDDO
ELSE
  ! No vector field provided
  IUVG = 0
  ISPUV = 0
  ZPGPUV=>NULL()
  ZPSPVOR=>NULL()
  ZPSPDIV=>NULL()
ENDIF

! 2. scalar fields transformation

! Preliminary checks

IF (PRESENT(YDFSPSCALAR) .NEQV. PRESENT(YDFSCALAR)) CALL ABOR1("[INV_TRANS_FIELD_API]  YDFSPSCALAR and YDFSCALAR must be provided together")

IF (PRESENT(YDFSPSCALAR)) THEN

  IF ((SIZE(YDFSPSCALAR)/= SIZE(YDFSCALAR))) CALL ABOR1("[INV_TRANS_FIELD_API] Inconsistent size for YDFSPSCALAR and YDFSCALAR")

  ! Convert list of spectral scalar fields of any domension into a list of 2d fields
  IFLDSPSC = LS_COUNT(YDFSPSCALAR)
  ALLOCATE(YLSPVSCALAR(IFLDSPSC))

  ALLOCATE(YLGVSCALAR(LG_COUNT(YDFSCALAR)))

  IFLDXG = SIZE(YLGVSCALAR) ! NUMBER OF OUTPUT SCALAR FIELDS IN GRID SPACE
  ! count the number of fields present on the processor
  C = LS(YLSPVSCALAR, YDFSPSCALAR, LDACC,.TRUE.)
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
  ALLOCATE(ZPSPSC2(IFLDXL,KSPEC))

! Allocate scalar field array in grid space
  ALLOCATE(ZPGP2(KPROMA,IFLDXG * ISCDIM,KGPBLKS))

! allocate 'b-set' for scalar fields
  ALLOCATE(IVSETSC2(IFLDXG))

  ! temporary copies on gpu
  if ( LDACC ) THEN
    !$ACC ENTER DATA CREATE(ZPSPSC2,ZPGP2)
  ENDIF

  ! Copy list of of spectral scalar fields into temporary arrays (1d copy thanks to field_view)

  ID = 1
  DO JFLD = 1,IFLDSPSC
    IF (ASSOCIATED(YLSPVSCALAR(JFLD)%P)) THEN
      ZZ1_1=>YLSPVSCALAR(JFLD)%P
      IF (LDACC) THEN
        !$ACC KERNELS PRESENT(ZPSPSC2,ZZ1_1)
        ZPSPSC2(ID,:) = ZZ1_1(:)
        !$ACC END KERNELS
      ELSE
        ZPSPSC2(ID,:) = ZZ1_1(:)
      ENDIF
      ID = ID + 1
    ENDIF
  ENDDO
  
  if (LDACC) THEN 
    !$ACC UPDATE SELF(ZPSPSC2)
  endif

 ! compute ´b-set´ for scalar-fields
  C = LG(YLGVSCALAR,YDFSCALAR, LDACC,.TRUE.)
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
    IF (PRESENT (FSPGL_PROC)) THEN
        CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                     & PSPSC2=ZPSPSC2,PGP2=ZPGP2,KVSETSC2=IVSETSC2, &
                     & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                     & KPROMA=KPROMA,FSPGL_PROC=FSPGL_PROC)
    ELSE
        CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                     & PSPSC2=ZPSPSC2,PGP2=ZPGP2, KVSETSC2=IVSETSC2, &
                     & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                     & KPROMA=KPROMA)
    ENDIF
ELSE IF (ASSOCIATED(ZPGP2)) THEN
    IF (PRESENT (FSPGL_PROC)) THEN
        CALL INV_TRANS(PSPSC2=ZPSPSC2,PGP2=ZPGP2,KVSETSC2=IVSETSC2, &
                     & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                     & KPROMA=KPROMA,FSPGL_PROC=FSPGL_PROC)
    ELSE
        CALL INV_TRANS(PSPSC2=ZPSPSC2,PGP2=ZPGP2, KVSETSC2=IVSETSC2, &
                     & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                     & KPROMA=KPROMA)
    ENDIF
ELSE IF (ASSOCIATED(ZPGPUV)) THEN
    IF (PRESENT (FSPGL_PROC)) THEN
        CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                     & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                     & KPROMA=KPROMA,FSPGL_PROC=FSPGL_PROC)
    ELSE
        CALL INV_TRANS(PSPVOR=ZPSPVOR,PSPDIV=ZPSPDIV,PGPUV=ZPGPUV,KVSETUV=IVSETUV, &
                     & LDSCDERS=LLSCDERS, LDVORGP=LLVORGP, LDDIVGP=LLDIVGP, LDUVDER=LLUVDER,  &
                     & KPROMA=KPROMA)
    ENDIF
ENDIF

! 4. Copy back temporary array data into grid-point fields

! remove garbage at the end of arrays
IEND = KGPTOT - KPROMA * (KGPBLKS - 1)
IF (LDACC) THEN
  IF (IUVG>0) THEN
    !$ACC UPDATE DEVICE(ZPGPUV)
    !$ACC KERNELS PRESENT(ZPGPUV)
    ZPGPUV (IEND+1:, :, :, KGPBLKS) = 0
    !$ACC END KERNELS
  ENDIF
  IF (IFLDXG>0) THEN
    !$ACC UPDATE DEVICE(ZPGP2)
    !$ACC KERNELS PRESENT(ZPGP2)
     ZPGP2 (IEND+1:, :, KGPBLKS) = 0
    !$ACC END KERNELS
  ENDIF
ELSE
  IF (IUVG>0) ZPGPUV (IEND+1:, :, :, KGPBLKS) = 0
  IF (IFLDXG>0)  ZPGP2 (IEND+1:, :, KGPBLKS) = 0
ENDIF

! copy vector fields

IF (IUVG>0) THEN

  IOFFSET = 0
  ! copy vorticity
  IF (LLVORGP) THEN
      C = LG(YLGVVOR,YDFVOR, LDACC, .FALSE.)
      DO JFLD=1,IUVG
        DO JLEV=1,KFLEVG
          ID = JLEV + (JFLD -1) * KFLEVG
          ZZ2_1=>YLGVVOR(ID)%P
          IF (LDACC) THEN
            !$ACC KERNELS PRESENT(ZPGPUV,ZZ2_1)
             ZZ2_1(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
            !$ACC END KERNELS
          ELSE
            ZZ2_1(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
          ENDIF
        ENDDO
      ENDDO

    IOFFSET = IOFFSET + 1
  ENDIF

  ! copy divergence
  IF (LLDIVGP) THEN

      C = LG(YLGVDIV,YDFDIV, LDACC, .FALSE.)
      DO JFLD=1,IUVG
        DO JLEV=1,KFLEVG
        ID = JLEV + (JFLD -1) * KFLEVG
        ZZ2_1=>YLGVDIV(ID)%P
        IF (LDACC) THEN
          !$ACC KERNELS PRESENT(ZPGPUV,ZZ2_1)
          ZZ2_1(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
          !$ACC END KERNELS
        ELSE
          ZZ2_1(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
        ENDIF
        ENDDO
      ENDDO

    IOFFSET = IOFFSET + 1
  ENDIF

  ! copy u and v
  C = LG(YLGVU, YDFU, LDACC, .FALSE.)
  C = LG(YLGVV, YDFV, LDACC, .FALSE.)


  DO JFLD=1,IUVG
    DO JLEV=1,KFLEVG
      ID = JLEV + (JFLD -1) * KFLEVG
      ZZ2_1=>YLGVU(ID)%P
      ZZ2_2=>YLGVV(ID)%P
      IF (LDACC) THEN
        !$ACC KERNELS PRESENT(ZPGPUV,ZZ2_1,ZZ2_2)
        ZZ2_1(:,:) =  ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:)
        ZZ2_2(:,:) =  ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:)
        !$ACC END KERNELS
      ELSE
        ZZ2_1(:,:) =  ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:)
        ZZ2_2(:,:) =  ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:)
      ENDIF
    ENDDO
  ENDDO

  IOFFSET = IOFFSET + 2

  ! copy u and v derivatives
  IF (LLUVDER) THEN
    C = LG(YLGVU_EW,YDFU_EW, LDACC, .FALSE.)
    C = LG(YLGVV_EW,YDFV_EW, LDACC, .FALSE.)

    DO JFLD=1,IUVG
      DO JLEV=1,KFLEVG
        ID = JLEV + (JFLD -1) * KFLEVG
        ZZ2_1=>YLGVU_EW(ID)%P
        ZZ2_2=>YLGVV_EW(ID)%P
        IF (LDACC) THEN
          !$ACC KERNELS PRESENT(ZPGPUV,ZZ2_1,ZZ2_2)
          ZZ2_1(:,:) =  ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:)
          ZZ2_2(:,:) =  ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:)
          !$ACC END KERNELS
        ELSE
          ZZ2_1(:,:) =  ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:)
          ZZ2_2(:,:) =  ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF

IF (IFLDXG > 0) THEN
  ! copy spectral scalar fields
    C = LG(YLGVSCALAR,YDFSCALAR, LDACC,.FALSE.)
    DO JFLD=1, IFLDXG
      ZZ2_1=>YLGVSCALAR(JFLD)%P(:,:)
      IF (LDACC) THEN
        !$ACC KERNELS PRESENT(ZPGP2,ZZ2_1)
        ZZ2_1(:,:) = ZPGP2(:,JFLD,:)
        !$ACC END KERNELS
      ELSE
        ZZ2_1(:,:) = ZPGP2(:,JFLD,:)
      ENDIF
    ENDDO

  ! copy spectral scalar fields derivatives

  IF (LLSCDERS) THEN
    C = LG(YLGVSCALAR_NS, YDFSCALAR_NS, LDACC, .FALSE.)
    C = LG(YLGVSCALAR_EW, YDFSCALAR_EW, LDACC, .FALSE.)

    DO JFLD=1,IFLDXG
        ZZ2_1=>YLGVSCALAR_NS(JFLD)%P
        ZZ2_2=>YLGVSCALAR_EW(JFLD)%P
        IF (LDACC) THEN
          !$ACC KERNELS PRESENT(ZPGP2,ZZ2_1,ZZ2_2)
          ZZ2_1(:,:) = ZPGP2(:, JFLD+IFLDXG,:)
          ZZ2_2(:,:) = ZPGP2(:, JFLD+(2*IFLDXG),:)
          !$ACC END KERNELS
        ELSE
          ZZ2_1(:,:) = ZPGP2(:, JFLD+IFLDXG,:)
          ZZ2_2(:,:) = ZPGP2(:, JFLD+(2*IFLDXG),:)
        ENDIF
      ENDDO

  ENDIF
ENDIF

! 5. Final cleanup

! delete temporary arrays
if ( LDACC ) THEN
  IF (ASSOCIATED(ZPSPVOR))  THEN
      !$ACC EXIT DATA DELETE(ZPSPVOR)
  ENDIF
  IF (ASSOCIATED(ZPSPDIV)) THEN
   !$ACC EXIT DATA DELETE(ZPSPDIV)
  ENDIF
  IF (ASSOCIATED(ZPGPUV)) THEN
     !$ACC EXIT DATA DELETE(ZPGPUV)
  ENDIF
  IF (ASSOCIATED(ZPSPSC2)) THEN
   !$ACC EXIT DATA DELETE(ZPSPSC2)
  ENDIF
  IF (ASSOCIATED(ZPGP2)) THEN
     !$ACC EXIT DATA DELETE(ZPGP2)
  ENDIF
ENDIF

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
