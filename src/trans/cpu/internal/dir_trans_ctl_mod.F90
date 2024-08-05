! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIR_TRANS_CTL_MOD
CONTAINS
SUBROUTINE DIR_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)

!**** *DIR_TRANS_CTL* - Control routine for direct spectral transform.

!     Purpose.
!     --------
!        Control routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS_CTL(...)

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

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

USE PARKIND1,          ONLY: JPIM, JPRB
USE TPM_GEN,           ONLY: NPROMATR, NOUT
USE SHUFFLE_MOD,       ONLY: SHUFFLE
USE FIELD_SPLIT_MOD,   ONLY: FIELD_SPLIT
USE LTDIR_CTL_MOD,     ONLY: LTDIR_CTL
USE FTDIR_CTL_MOD,     ONLY: FTDIR_CTL_COMM_SEND, FTDIR_CTL_COMP
USE ABORT_TRANS_MOD,   ONLY: ABORT_TRANS
USE OVERLAP_TYPES_MOD, ONLY: BATCH

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

! Local variables
INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)

TYPE(BATCH) :: NEW_BATCH
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV
INTEGER(KIND=JPIM) :: NDONE

!     ------------------------------------------------------------------

! Perform transform

WRITE(NOUT,*) "KF_UV_G = ", KF_UV_G
WRITE(NOUT,*) "KF_SCALARS_G = ", KF_SCALARS_G
WRITE(NOUT,*) "KF_GP = ", KF_GP
WRITE(NOUT,*) "KF_FS = ", KF_FS
WRITE(NOUT,*) "KF_UV = ", KF_UV
WRITE(NOUT,*) "KF_SCALARS = ", KF_SCALARS

IF (NPROMATR > 0 .AND. KF_GP > NPROMATR) THEN

  CALL SHUFFLE(KF_UV_G, KF_SCALARS_G, ISHFUV_G, IVSETUV, ISHFSC_G, IVSETSC, KVSETUV, KVSETSC)

  IBLKS=(KF_GP-1)/NPROMATR+1

  JBLK = 1
  NDONE = 0

  DO WHILE (NDONE < IBLKS)
    NEW_BATCH = BATCH(JBLK, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC)

    CALL FIELD_SPLIT(JBLK, KF_GP, KF_UV_G, IVSETUV, IVSETSC, ISTUV_G, IENUV_G, IF_UV_G, ISTSC_G, &
      &              IENSC_G, IF_SCALARS_G, ISTUV, IENUV, IF_UV, ISTSC, IENSC, IF_SCALARS)

    IF_FS = 2*IF_UV + IF_SCALARS
    IF_GP = 2*IF_UV_G+IF_SCALARS_G

    DO JFLD=1,IF_UV_G
      IPTRGP(JFLD) = ISHFUV_G(ISTUV_G+JFLD-1)
      IPTRGP(JFLD+IF_UV_G) = KF_UV_G+ISHFUV_G(ISTUV_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_SCALARS_G
      IPTRGP(JFLD+2*IF_UV_G) = 2*KF_UV_G+ISHFSC_G(ISTSC_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_UV
      IPTRSPUV(JFLD) = ISTUV+JFLD-1
    ENDDO
    DO JFLD=1,IF_SCALARS
      IPTRSPSC(JFLD) = ISTSC+JFLD-1
    ENDDO

    IF (NEW_BATCH%NF_UV_G > 0 .AND. NEW_BATCH%NF_SCALARS_G > 0) THEN
      CALL FTDIR_CTL_COMM_SEND(NEW_BATCH%NF_UV_G, NEW_BATCH%NF_SCALARS_G, NEW_BATCH%NF_GP, NEW_BATCH%NF_FS, &
        &                      KVSETUV=IVSETUV(ISTUV_G:IENUV_G), KVSETSC=IVSETSC(ISTSC_G:IENSC_G), &
        &                      KPTRGP=IPTRGP, PGP=PGP)
    ELSEIF (NEW_BATCH%NF_UV_G > 0) THEN
      CALL FTDIR_CTL_COMM_SEND(NEW_BATCH%NF_UV_G, NEW_BATCH%NF_SCALARS_G, NEW_BATCH%NF_GP, NEW_BATCH%NF_FS, &
        &                      KVSETUV=IVSETUV(ISTUV_G:IENUV_G), KPTRGP=IPTRGP, PGP=PGP)
    ELSEIF (NEW_BATCH%NF_SCALARS_G > 0) THEN
      CALL FTDIR_CTL_COMM_SEND(NEW_BATCH%NF_UV_G, NEW_BATCH%NF_SCALARS_G, NEW_BATCH%NF_GP, NEW_BATCH%NF_FS, &
        &                      KVSETSC=IVSETSC(ISTSC_G:IENSC_G), KPTRGP=IPTRGP, PGP=PGP)
    ENDIF
    CALL FTDIR_CTL_COMP(NEW_BATCH%NF_FS)
    CALL LTDIR_CTL(NEW_BATCH%NF_FS, NEW_BATCH%NF_UV, NEW_BATCH%NF_SCALARS, PSPVOR=PSPVOR, PSPDIV=PSPDIV, PSPSCALAR=PSPSCALAR, &
      &            KFLDPTRUV=IPTRSPUV, KFLDPTRSC=IPTRSPSC)

    !CALL NEW_BATCH%START_COMM(PGP, PSPVOR, PSPDIV, PSPSCALAR)
    JBLK = JBLK + 1
    NDONE = NDONE + 1
  ENDDO
ELSE

  ! No splitting of fields, transform done in one go
  CALL ABORT_TRANS("DIR_TRANS_CTL: NPROMATR = 0 feature disabled for overlap version")

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_CTL
END MODULE DIR_TRANS_CTL_MOD
