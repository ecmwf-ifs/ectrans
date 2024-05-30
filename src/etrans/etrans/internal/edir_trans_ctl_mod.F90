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
!        G. Radnoti 01-03-13 adaptation to aladin
!     01-08-28 : G. Radnoti & R. El Khatib Fix for NPROMATR /= 0
!     02-09-30 : P. Smolikova AUX_PROC for d4 in NH
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!    ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : NPROMATR
!USE TPM_TRANS
!USE TPM_DISTR

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
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

INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_GPB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! Perform transform

IF (LHOOK) CALL DR_HOOK('EDIR_TRANS_CTL_MOD:EDIR_TRANS_CTL',0,ZHOOK_HANDLE)
IF_GPB = 2*KF_UV_G+KF_SCALARS_G
IF(NPROMATR > 0 .AND. IF_GPB > NPROMATR) THEN

  ! Fields to be split into packets

  CALL SHUFFLE(KF_UV_G,KF_SCALARS_G,ISHFUV_G,IVSETUV,ISHFSC_G,IVSETSC, &
   & KVSETUV,KVSETSC)

  IBLKS=(IF_GPB-1)/NPROMATR+1

  DO JBLK=1,IBLKS

    CALL FIELD_SPLIT(JBLK,KF_GP,KF_UV_G,IVSETUV,IVSETSC,&
     & ISTUV_G,IENUV_G,IF_UV_G,ISTSC_G,IENSC_G,IF_SCALARS_G,&
     & ISTUV,IENUV,IF_UV,ISTSC,IENSC,IF_SCALARS)

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

    IF(IF_UV_G > 0 .AND. IF_SCALARS_G > 0) THEN
      CALL EFTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,IF_GPB,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,PGP=PGP)
    ELSEIF(IF_UV_G > 0) THEN
      CALL EFTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,IF_GPB,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KPTRGP=IPTRGP,PGP=PGP)
    ELSEIF(IF_SCALARS_G > 0) THEN
      CALL EFTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,IF_GPB,&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,PGP=PGP,&
       & AUX_PROC=AUX_PROC)
    ENDIF
    CALL ELTDIR_CTL(IF_FS,IF_UV,IF_SCALARS, &
     & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
     & KFLDPTRUV=IPTRSPUV,KFLDPTRSC=IPTRSPSC,&
     & PSPMEANU=PMEANU,PSPMEANV=PMEANV,AUX_PROC=AUX_PROC)
  ENDDO
ELSE

  ! No splitting of fields, transform done in one go

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

ENDIF
IF (LHOOK) CALL DR_HOOK('EDIR_TRANS_CTL_MOD:EDIR_TRANS_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EDIR_TRANS_CTL
END MODULE EDIR_TRANS_CTL_MOD
