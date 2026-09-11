! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIR_TRANS_VIEW_CTL_MOD
CONTAINS
SUBROUTINE DIR_TRANS_VIEW_CTL(KPROMA,KGPBLKS, &
                            & YDSPVVOR, YDSPVDIV, YDSPVSCALAR, &
                            & YDGVU,YDGVV,YDGVSCALAR)

!**** *DIR_TRANS_VIEW_CTL* - Control routine for direct spectral transform.

!     Purpose.
!     --------
!        Control routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS_VIEW_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity
!     PSPDIV(:,:)  - spectral divergence
!     PSPSCALAR(:,:) - spectral scalarvalued fields
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
!     PGP(:,:,:)  - gridpoint fields

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTDIR_CTL   - control of Legendre transform
!                 FTDIR_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NPROMATR
!USE TPM_TRANS
!USE TPM_DISTR

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
USE LTDIR_VIEW_CTL_MOD   ,ONLY : LTDIR_VIEW_CTL
USE FTDIR_VIEW_CTL_MOD   ,ONLY : FTDIR_VIEW_CTL
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW!
USE TPM_TRANS       ,ONLY : NPROMA, NGPBLKS
USE TPM_DISTR       ,ONLY : MYSETV

IMPLICIT NONE

! Declaration of arguments
INTEGER(KIND=JPIM) :: KPROMA, KGPBLKS
TYPE(SPEC_VIEW) :: YDSPVVOR(:), YDSPVDIV(:)
TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)

TYPE(GRID_VIEW) :: YDGVU(:),YDGVV(:)
TYPE(GRID_VIEW) :: YDGVSCALAR(:)

! Local variables

INTEGER(KIND=JPIM) :: IF_FS,IF_GP, IF_UV, J
TYPE(GRID_VIEW),ALLOCATABLE :: YLGP(:)
!     ------------------------------------------------------------------

! Perform transform
NPROMA = KPROMA
NGPBLKS = KGPBLKS

IF(NPROMATR > 0) THEN
  CALL ABOR1('DIR_TRANS_VIEW_CTL NPROMATR not implemented')
ELSE

  YLGP = [YDGVU, YDGVV,YDGVSCALAR]
  IF_GP = SIZE(YLGP)
  IF_FS = 0
  DO J = 1, IF_GP
    IF (YLGP(J)%IVSET == MYSETV .OR. YLGP(J)%IVSET == -1) IF_FS = IF_FS + 1
  ENDDO

  CALL FTDIR_VIEW_CTL(IF_GP,IF_FS, YDGP=YLGP)

  CALL LTDIR_VIEW_CTL(IF_FS,YDSPVVOR=YDSPVVOR,YDSPVDIV=YDSPVDIV,YDSPVSCALAR=YDSPVSCALAR)

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_VIEW_CTL
END MODULE DIR_TRANS_VIEW_CTL_MOD
