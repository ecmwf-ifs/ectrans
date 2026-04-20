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

TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGU(:), YLFVGV(:)
TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGVOR(:), YLFVGDIV(:)
TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGU_EW(:), YLFVGV_EW(:)

TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGSCALAR(:)

TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGSCALAR_NS(:)
TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVGSCALAR_EW(:)

TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVSVOR(:), YLFVSDIV(:)
TYPE(FIELD_VIEW), ALLOCATABLE  :: YLFVSSCALAR(:)

INTEGER, PARAMETER :: FIELD_TYPE = 0

INTEGER(KIND=JPIM) :: IFIELD
INTEGER(KIND=JPIM) :: NFIELD_UV, NFIELD_SCALAR
INTEGER(KIND=JPIM) :: NFIELD_UVDER, NFIELD_DIVGP, NFIELD_SCDER, NFIELD_VORGP

REAL(KIND=JPHOOK)           :: ZHOOK_HANDLE

#include "inv_trans.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_API',0,ZHOOK_HANDLE)


! -------------------------------------------------------------------------------------------------------------
! Convert FIELD_API input fields into FIELD_VIEW, an intermediate representation of fields specific to ectrans
NFIELD_UV = 0
NFIELD_SCALAR = 0
NFIELD_UVDER = 0
NFIELD_DIVGP = 0
NFIELD_SCDER = 0
NFIELD_VORGP = 0

IF (PRESENT(YDFU) .AND. PRESENT(YDFV)) THEN
  NFIELD_UV = SIZE(YDFU)
  ENDIF
IF (PRESENT(YDFSCALAR)) THEN
  NFIELD_SCALAR = SIZE(YDFSCALAR)
ENDIF
  IF (PRESENT(YDFU_EW) .AND. PRESENT(YDFV_EW))    THEN
  NFIELD_UVDER = NFIELD_UV
 ENDIF

  IF (PRESENT(YDFDIV)) THEN
    NFIELD_DIVGP = NFIELD_UV
  ENDIF

  IF (PRESENT(YDFVOR)) THEN
    NFIELD_VORGP = NFIELD_UV
  ENDIF

  ! allocate temporary vector field arrays in spectral space
  ALLOCATE(ZPSPVOR(IFLDSPVOR,NSPEC2))
  ALLOCATE(ZPSPDIV(IFLDSPVOR,NSPEC2))

  ! allocate temporary vector field array in grid space
  ALLOCATE(ZPGPUV(NPROMA,NFLEVG, IUVG * IUVDIM,NBLK))

  ! allocate 'b-set' for vector fields
  ALLOCATE(IVSETUV(NFLEVG))

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
  C = LG(YDFU, YLGVU, LDACC, .TRUE.)
  DO JFLD=1,IUVG
    DO JLEV=1,NFLEVG
     ID = JLEV + (JFLD -1) * NFLEVG
     IF (JFLD .EQ. 1) IVSETUV(JLEV) = YLGVU(ID)%IVSET
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
    NFIELD_SCDER = NFIELD_SCALAR
ENDIF

ALLOCATE(YLFVSVOR(NFIELD_UV))
ALLOCATE(YLFVSDIV(NFIELD_UV))
ALLOCATE(YLFVSSCALAR(NFIELD_SCALAR))

ALLOCATE(YLFVGU(NFIELD_UV))
ALLOCATE(YLFVGV(NFIELD_UV))
ALLOCATE(YLFVGVOR(NFIELD_VORGP))
ALLOCATE(YLFVGDIV(NFIELD_DIVGP))
ALLOCATE(YLFVGSCALAR(NFIELD_SCALAR))
ALLOCATE(YLFVGU_EW(NFIELD_UVDER))
ALLOCATE(YLFVGV_EW(NFIELD_UVDER))
ALLOCATE(YLFVGSCALAR_NS(NFIELD_SCDER))
ALLOCATE(YLFVGSCALAR_EW(NFIELD_SCDER))

DO IFIELD=1,NFIELD_UV
  CALL MAKE_FIELD_VIEW(YLFVGU(IFIELD),   YDFU(IFIELD),     FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVGV(IFIELD),   YDFV(IFIELD),     FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVSVOR(IFIELD), YDFSPVOR(IFIELD), FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVSDIV(IFIELD), YDFSPDIV(IFIELD), FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
ENDDO

DO IFIELD=1,NFIELD_UVDER
   CALL MAKE_FIELD_VIEW(YLFVGU_EW(IFIELD),   YDFU_EW(IFIELD),     FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
   CALL MAKE_FIELD_VIEW(YLFVGV_EW(IFIELD),   YDFV_EW(IFIELD),     FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
ENDDO

DO IFIELD=1,NFIELD_DIVGP
   CALL MAKE_FIELD_VIEW(YLFVGDIV(IFIELD),   YDFDIV(IFIELD),     FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
ENDDO

DO IFIELD=1,NFIELD_VORGP
  CALL MAKE_FIELD_VIEW(YLFVGVOR(IFIELD),   YDFVOR(IFIELD),     FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
ENDDO

DO IFIELD=1,NFIELD_SCALAR
  CALL MAKE_FIELD_VIEW(YLFVGSCALAR(IFIELD), YDFSCALAR(IFIELD),   FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVSSCALAR(IFIELD), YDFSPSCALAR(IFIELD), FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
ENDDO

DO IFIELD=1,NFIELD_SCDER
  CALL MAKE_FIELD_VIEW(YLFVGSCALAR_NS(IFIELD), YDFSCALAR_NS(IFIELD), FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
  CALL MAKE_FIELD_VIEW(YLFVGSCALAR_EW(IFIELD), YDFSCALAR_EW(IFIELD), FIELD_TYPE, LDACC=.FALSE., LDRDONLY=.TRUE.)
ENDDO

CALL INV_TRANS_FIELD_VIEW(KRESOL, &
                        & YLFVSVOR,YLFVSDIV,YLFVSSCALAR, &
                        & YLFVGU, YLFVGV, YLFVGVOR,YLFVGDIV,YLFVGSCALAR, &
                        & YLFVGU_EW, YLFVGV_EW, YLFVGSCALAR_NS, YLFVGSCALAR_EW,&
                        & KGPTOT, &
                        & FSPGL_PROC)

DEALLOCATE(YLFVSVOR)
DEALLOCATE(YLFVSDIV)
DEALLOCATE(YLFVSSCALAR)

DEALLOCATE(YLFVGU)
DEALLOCATE(YLFVGV)
DEALLOCATE(YLFVGSCALAR)
DEALLOCATE(YLFVGVOR)
DEALLOCATE(YLFVGDIV)
DEALLOCATE(YLFVGU_EW)
DEALLOCATE(YLFVGV_EW)
DEALLOCATE(YLFVGSCALAR_NS)
DEALLOCATE(YLFVGSCALAR_EW)

IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_API',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_FIELD_API
