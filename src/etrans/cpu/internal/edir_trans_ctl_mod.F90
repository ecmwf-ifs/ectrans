! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EDIR_TRANS_CTL_MOD
CONTAINS
SUBROUTINE EDIR_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2,&
 & PMEANU,PMEANV,AUX_PROC)

!**** *EDIR_TRANS_CTL* - Control routine for direct spectral transform.

!     Purpose.
!     --------
!        Control routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL EDIR_TRANS_CTL(...)

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
!    PMEANU,PMEANV - mean winds
!    AUX_PROC        - optional external procedure for biperiodization of
!            aux.fields
!     PGP(:,:,:)  - gridpoint fields

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):

!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  LTDIR_CTL   - control of Legendre transform
!     ----------  FTDIR_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03
!        G. Radnoti 01-03-13 adaptation to aladin
!     01-08-28 : G. Radnoti & R. El Khatib Fix for NPROMATR /= 0
!     02-09-30 : P. Smolikova AUX_PROC for d4 in NH
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!    ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE ELTDIR_CTL_MOD  ,ONLY : ELTDIR_CTL
USE EFTDIR_CTL_MOD  ,ONLY : EFTDIR_CTL

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
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
EXTERNAL AUX_PROC
OPTIONAL AUX_PROC

! Local variables
INTEGER(KIND=JPIM) :: IF_GPB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! Perform transform

IF (LHOOK) CALL DR_HOOK('EDIR_TRANS_CTL_MOD:EDIR_TRANS_CTL',0,ZHOOK_HANDLE)

IF_GPB = 2*KF_UV_G+KF_SCALARS_G
CALL EFTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,IF_GPB,&
  & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
  & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
  & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2,&
  & AUX_PROC=AUX_PROC)

CALL ELTDIR_CTL(KF_FS,KF_UV,KF_SCALARS, &
  & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
  & PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2,&
  & PSPMEANU=PMEANU,PSPMEANV=PMEANV,&
  & AUX_PROC=AUX_PROC)

IF (LHOOK) CALL DR_HOOK('EDIR_TRANS_CTL_MOD:EDIR_TRANS_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EDIR_TRANS_CTL
END MODULE EDIR_TRANS_CTL_MOD
