! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EINV_TRANS_CTLAD_MOD
CONTAINS
SUBROUTINE EINV_TRANS_CTLAD(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT,&
 & KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2,&
 & PMEANU,PMEANV)

!**** *EINV_TRANS_CTLAD* - Control routine for inverse spectral transform adj.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL EINV_TRANS_CTLAD(...)

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

!     Externals.  LTINV_CTLAD   - control of Legendre transform
!     ----------  FTINV_CTLAD   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE ELTINV_CTLAD_MOD ,ONLY : ELTINV_CTLAD
USE EFTINV_CTLAD_MOD ,ONLY : EFTINV_CTLAD
!

IMPLICIT NONE

! Declaration of arguments
!
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC2(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPUV(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3A(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3B(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP2(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PMEANU(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PMEANV(:)

! Local variables
INTEGER(KIND=JPIM) :: IF_GPB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! Perform transform

IF (LHOOK) CALL DR_HOOK('EINV_TRANS_CTLAD_MOD:EINV_TRANS_CTLAD',0,ZHOOK_HANDLE)

IF_GPB = 2*KF_UV_G+KF_SCALARS_G
CALL EFTINV_CTLAD(KF_UV_G,KF_SCALARS_G,&
  & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,&
  & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
  & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
  & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)

CALL ELTINV_CTLAD(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS, &
  & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
  & PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2,&
  & PSPMEANU=PMEANU,PSPMEANV=PMEANV )

IF (LHOOK) CALL DR_HOOK('EINV_TRANS_CTLAD_MOD:EINV_TRANS_CTLAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EINV_TRANS_CTLAD
END MODULE EINV_TRANS_CTLAD_MOD
