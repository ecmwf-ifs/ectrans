! (C) Copyright 2001- ECMWF.
! (C) Copyright 2022- NVIDIA.
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
USE PARKIND_ECTRANS  ,ONLY : JPRBT

USE TPM_GEN         ,ONLY : NPROMATR, NOUT
USE TPM_TRANS       ,ONLY : NF_SC2, NF_SC3A, NF_SC3B
USE TPM_DISTR, ONLY: NPROC

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
USE LTDIR_CTL_MOD   ,ONLY : LTDIR_CTL
USE FTDIR_CTL_MOD   ,ONLY : FTDIR_CTL
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE TRLTOM_MOD      ,ONLY : TRLTOM, TRLTOM_CUDAAWARE
!

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
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_GPB

REAL(KIND=JPRBT), ALLOCATABLE :: FOUBUF_IN(:)
REAL(KIND=JPRBT), ALLOCATABLE :: FOUBUF(:)


!     ------------------------------------------------------------------

! Perform transform

IF_GPB = 2*KF_UV_G+KF_SCALARS_G
IF(NPROMATR > 0 .AND. IF_GPB > NPROMATR) THEN
    PRINT *, "ERROR, not implemented right now"
    stop 4
ELSE

  ! No splitting of fields, transform done in one go
  CALL FTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,FOUBUF_IN,&
   & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
   & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
   & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)

  CALL GSTATS(153,0)
  IF (NPROC > 1) THEN
#ifdef USE_CUDA_AWARE_MPI_FT
    WRITE(NOUT,*) 'ltdir_ctl:TRLTOM_CUDAAWARE'
    CALL TRLTOM_CUDAAWARE(FOUBUF_IN,FOUBUF,KF_FS)
#else
    CALL TRLTOM(FOUBUF_IN,FOUBUF,KF_FS)
#endif
  ENDIF
  CALL GSTATS(153,1)

  CALL LTDIR_CTL(KF_FS,KF_UV,KF_SCALARS,FOUBUF, &
   & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
   & PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_CTL
END MODULE DIR_TRANS_CTL_MOD
