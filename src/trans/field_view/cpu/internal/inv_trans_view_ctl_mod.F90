! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS_VIEW_CTL_MOD
CONTAINS
SUBROUTINE INV_TRANS_VIEW_CTL(KPROMA,KGPBLKS, &
                            & YDSPVVOR, YDSPVDIV, YDSPVSCALAR, &
                            & YDGVU,YDGVV,&
                            & YDGVVOR,YDGVDIV,&
                            & YDGVSCALAR,&
                            & YDGVU_EW,YDGVV_EW,&
                            & YDGVSCALAR_NS, YDGVSCALAR_EW,&
                            & FSPGL_PROC)

!**** *INV_TRANS_VIEW_CTL* - Control routine for inverse spectral transform.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_VIEW_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_OUT_LT    - total number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition
!     PGP(:,:,:)  - gridpoint fields (output)

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):

!       vorticity     : KF_UV_G fields
!       divergence    : KF_UV_G fields
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields
!       N-S derivative of scalar fields : KF_SCALARS_G fields
!       E-W derivative of u : KF_UV_G fields
!       E-W derivative of v : KF_UV_G fields
!       E-W derivative of scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTINV_CTL   - control of Legendre transform
!                 FTINV_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NPROMATR
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP, NPROMA, NGPBLKS
USE TPM_DISTR       ,ONLY : MYSETV

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
USE LTINV_VIEW_CTL_MOD   ,ONLY : LTINV_VIEW_CTL
USE FTINV_VIEW_CTL_MOD   ,ONLY : FTINV_VIEW_CTL
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW

IMPLICIT NONE

#include "fspgl_intf.h"

! Declaration of arguments
INTEGER(KIND=JPIM) :: KPROMA, KGPBLKS
TYPE(SPEC_VIEW) :: YDSPVVOR(:), YDSPVDIV(:)
TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)

TYPE(GRID_VIEW) :: YDGVU(:),YDGVV(:)
TYPE(GRID_VIEW) :: YDGVVOR(:),YDGVDIV(:)
TYPE(GRID_VIEW) :: YDGVSCALAR(:)

TYPE(GRID_VIEW) :: YDGVU_EW(:),YDGVV_EW(:)
TYPE(GRID_VIEW) :: YDGVSCALAR_NS(:), YDGVSCALAR_EW(:)

PROCEDURE (FSPGL_INTF), POINTER, OPTIONAL, INTENT(IN)  :: FSPGL_PROC

INTEGER(KIND=JPIM) :: IF_UV,IF_UV_G,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP,IF_OUT_LT
INTEGER(KIND=JPIM) :: IF_SCDERS,IF_SCDERS_G,IF_UV_PAR
INTEGER(KIND=JPIM) :: IF_SC2_G,IF_SC3A_G2,IF_SC3A_G3,IF_SC3B_G2,IF_SC3B_G3
INTEGER(KIND=JPIM) :: J

NPROMA = KPROMA
NGPBLKS = KGPBLKS

IF_UV_G = SIZE(YDGVU)
IF_SCALARS_G = SIZE(YDGVSCALAR)
IF_SCDERS_G = SIZE(YDGVSCALAR_NS)

IF_UV = 0
DO J = 1, IF_UV_G
  IF (YDGVU(J)%IVSET == MYSETV .OR. YDGVU(J)%IVSET == -1) IF_UV = IF_UV + 1
ENDDO
IF_SCALARS = 0
DO J = 1, IF_SCALARS_G
  IF (YDGVSCALAR(J)%IVSET == MYSETV .OR. YDGVSCALAR(J)%IVSET == -1) IF_SCALARS = IF_SCALARS + 1
ENDDO
IF_SCDERS = 0
DO J = 1, IF_SCDERS_G
  IF (YDGVSCALAR_NS(J)%IVSET == MYSETV .OR. YDGVSCALAR_NS(J)%IVSET == -1) IF_SCDERS = IF_SCDERS + 1
ENDDO

LVORGP = SIZE(YDGVVOR) > 0
LDIVGP = SIZE(YDGVDIV) > 0
LUVDER = SIZE(YDGVU_EW) > 0 .OR. SIZE(YDGVV_EW) > 0
LSCDERS = SIZE(YDGVSCALAR_NS) > 0 .OR. SIZE(YDGVSCALAR_EW) > 0

IF_GP = SIZE(YDGVU) + SIZE(YDGVV) + SIZE(YDGVSCALAR) + SIZE(YDGVVOR) + SIZE(YDGVDIV) +&
        &SIZE(YDGVU_EW) + SIZE(YDGVV_EW) + SIZE(YDGVSCALAR_NS) + SIZE(YDGVSCALAR_EW)


!local number of fields coming out from inverse LT (per V-set)
IF_OUT_LT = 2*IF_UV + IF_SCALARS+IF_SCDERS

IF(IF_UV > 0 .AND. LVORGP) THEN
  IF_OUT_LT = IF_OUT_LT+IF_UV
ENDIF
IF(IF_UV > 0 .AND. LDIVGP) THEN
  IF_OUT_LT = IF_OUT_LT+IF_UV
ENDIF
IF_FS = IF_OUT_LT+IF_SCDERS
IF(IF_UV > 0 .AND. LUVDER) THEN
  IF_FS = IF_FS+2*IF_UV
ENDIF

!     ------------------------------------------------------------------

! Perform transform

IF(NPROMATR > 0) THEN
  CALL ABOR1('INV_TRANS_VIEW NPROMATR not implemented')
ELSE

  CALL LTINV_VIEW_CTL(IF_OUT_LT,&
   &YDSPVVOR=YDSPVVOR,YDSPVDIV=YDSPVDIV,YDSPVSCALAR=YDSPVSCALAR,&
   &FSPGL_PROC=FSPGL_PROC)

  CALL FTINV_VIEW_CTL(IF_GP, IF_FS, IF_OUT_LT, &
   & YDGVU=YDGVU,YDGVV=YDGVV,YDGVSCALAR=YDGVSCALAR,&
   & YDGVVOR=YDGVVOR,YDGVDIV=YDGVDIV,&
   & YDGVU_EW=YDGVU_EW,YDGVV_EW=YDGVV_EW, &
   & YDGVSCALAR_NS=YDGVSCALAR_NS,YDGVSCALAR_EW=YDGVSCALAR_EW)
ENDIF

!     ------------------------------------------------------------------
END SUBROUTINE INV_TRANS_VIEW_CTL
END MODULE INV_TRANS_VIEW_CTL_MOD
