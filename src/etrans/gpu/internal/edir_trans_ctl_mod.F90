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
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP, GROWING_ALLOCATION
!USE TPM_TRANS
USE TPM_DISTR, only : D, MYSETW


USE ELTDIR_MOD
USE TRLTOM_PACK_UNPACK, ONLY : TRLTOM_PACK, TRLTOM_PACK_HANDLE, PREPARE_TRLTOM_PACK
USE TRLTOM_MOD
USE EFTDIR_MOD
USE TRGTOL_MOD
USE BUFFERED_ALLOCATOR_MOD

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

REAL(KIND=JPRB), POINTER :: FOUBUF(:), FOUBUF_IN(:), PREEL(:)
INTEGER(KIND=JPIM) :: ILEI2, IDIM1
TYPE(BUFFERED_ALLOCATOR) :: ALLOCATOR
TYPE(ELTDIR_HANDLE) :: HELTDIR
TYPE(TRLTOM_HANDLE) :: HTRLTOM
TYPE(TRLTOM_PACK_HANDLE) :: HTRLTOM_PACK
TYPE(TRGTOL_HANDLE) :: HTRGTOL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!write (6,*) __FILE__, __LINE__; call flush(6)
!write (6,*) 'entered edir_trans_ctl'
!call flush(6)

! Perform transform

IF (LHOOK) CALL DR_HOOK('EDIR_TRANS_CTL_MOD:EDIR_TRANS_CTL',0,ZHOOK_HANDLE)
IF_GPB = 2*KF_UV_G+KF_SCALARS_G
IF(NPROMATR > 0) THEN
  print *, "This is currently not supported and/or tested (NPROMATR > 0)"
  stop 24
ENDIF


! Prepare everything
ALLOCATOR = MAKE_BUFFERED_ALLOCATOR()
HTRGTOL = PREPARE_TRGTOL(ALLOCATOR,KF_GP,KF_FS)  ! ZCOMBUFR, ZCOMBUFS and PREEL
HTRLTOM_PACK = PREPARE_TRLTOM_PACK(ALLOCATOR, KF_FS) ! FOUBUF_IN
HTRLTOM = PREPARE_TRLTOM(ALLOCATOR, KF_FS) ! FOUBUF
HELTDIR = PREPARE_ELTDIR(ALLOCATOR, KF_FS, KF_UV)

CALL INSTANTIATE_ALLOCATOR(ALLOCATOR, GROWING_ALLOCATION)

!write (6,*) __FILE__, __LINE__; call flush(6)

! from the PGP arrays to PREEL_REAL
CALL TRGTOL(ALLOCATOR,HTRGTOL,PREEL,KF_FS,KF_GP,KF_UV_G,KF_SCALARS_G,&
 & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
 & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
 & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)


!write (6,*) __FILE__, __LINE__; call flush(6)

IF (KF_FS > 0) THEN

!write (6,*) __FILE__, __LINE__; call flush(6)

  ! fourier transform from PREEL_REAL to PREEL_COMPLEX (in-place!)
  CALL GSTATS(1640,0)
  CALL EFTDIR(ALLOCATOR,PREEL,KF_FS,AUX_PROC=AUX_PROC)
  CALL GSTATS(1640,1)

!write (6,*) __FILE__, __LINE__; call flush(6)

  CALL GSTATS(153,0)
  CALL TRLTOM_PACK(ALLOCATOR,HTRLTOM_PACK,PREEL,FOUBUF_IN,KF_FS)    ! formerly known as efourier_out
!write (6,*) __FILE__, __LINE__; call flush(6)
  CALL TRLTOM(ALLOCATOR,HTRLTOM,FOUBUF_IN,FOUBUF,KF_FS)
  CALL GSTATS(153,1)

!write (6,*) __FILE__, __LINE__; call flush(6)

  CALL ELTDIR(ALLOCATOR,HELTDIR,KF_FS,KF_UV,KF_SCALARS,FOUBUF, &
        & PSPVOR,PSPDIV,PSPSCALAR,&
        & PSPSC3A,PSPSC3B,PSPSC2, &
        & PSPMEANU=PMEANU,PSPMEANV=PMEANV)

!write (6,*) 'leaving edir_trans_ctl'
!call flush(6)
!write (6,*) __FILE__, __LINE__; call flush(6)

ENDIF

IF (LHOOK) CALL DR_HOOK('EDIR_TRANS_CTL_MOD:EDIR_TRANS_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EDIR_TRANS_CTL
END MODULE EDIR_TRANS_CTL_MOD
