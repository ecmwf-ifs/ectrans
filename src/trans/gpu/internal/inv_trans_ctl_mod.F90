! (C) Copyright 2001- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS_CTL_MOD
CONTAINS
SUBROUTINE INV_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT,&
 & KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,FSPGL_PROC,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)

!**** *INV_TRANS_CTL* - Control routine for inverse spectral transform.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTL(...)

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


USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,   JPRBT

USE TPM_GEN         ,ONLY : NPROMATR, NOUT, NERR
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP, REUSE_PTR
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE TRMTOL_MOD
USE LTINV_MOD
USE FOURIER_IN_MOD
USE FSC_MOD
USE FTINV_MOD
USE TRLTOG_MOD
USE TPM_DISTR       ,ONLY : D
USE ALLOCATOR_MOD
!

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
REAL(KIND=JPRB) ,OPTIONAL   ,INTENT(OUT) :: PGP(:,:,:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP2(:,:,:)

! Local variables

REAL(KIND=JPRB), POINTER :: FOUBUF(:)
REAL(KIND=JPRB), POINTER :: FOUBUF_IN(:)

REAL(KIND=JPRBT), POINTER :: PREEL_REAL(:)
REAL(KIND=JPRBT), POINTER :: PREEL_COMPLEX(:)
INTEGER(KIND=JPIM) :: KUV_OFFSET, KSCALARS_OFFSET, KSCALARS_NSDER_OFFSET, &
    & KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET
INTEGER(KIND=JPIM) :: IF_LEG, IF_FOURIER

INTEGER(KIND=JPIM) :: IST
INTEGER(KIND=JPIM) :: IVSETUV(KF_UV_G)
INTEGER(KIND=JPIM) :: IVSETSC(KF_SCALARS_G)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: J3,JGL,IGL,IOFF,IFGP2,IFGP3A,IFGP3B,IGP3APAR,IGP3BPAR
INTEGER(KIND=JPIM) :: IFIRST

TYPE(BUFFERED_ALLOCATOR) :: ALLOCATOR
TYPE(LTINV_HANDLE) :: HLTINV
TYPE(TRMTOL_HANDLE) :: HTRMTOL
TYPE(FOURIER_IN_HANDLE) :: HFOURIER_IN
TYPE(FSC_HANDLE) :: HFSC
TYPE(FTINV_HANDLE) :: HFTINV
TYPE(TRLTOG_HANDLE) :: HTRLTOG

INTEGER(KIND=C_INT8_T), POINTER :: PTR(:)

!     ------------------------------------------------------------------

  IF(NPROMATR > 0) THEN
    print *, "This is currently not supported and/or tested (NPROMATR > 0j"
    stop 24
  ENDIF

  ! Compute Vertical domain decomposition

  ! Initialize potentially unset offsets
  KSCALARS_NSDER_OFFSET = -1
  KUV_EWDER_OFFSET = -1
  KSCALARS_EWDER_OFFSET = -1

  ! (note in ltinv we will initially start with a slightly different domain decomposition
  ! which always has vorticity and divergence because this is the actual input)
  IFIRST = 0
  IF (LVORGP) IFIRST = IFIRST + KF_UV ! Vorticity
  IF (LDIVGP) IFIRST = IFIRST + KF_UV ! Divergence
  KUV_OFFSET = IFIRST
  IFIRST = IFIRST + KF_UV ! U
  IFIRST = IFIRST + KF_UV ! V
  KSCALARS_OFFSET = IFIRST
  IFIRST = IFIRST + KF_SCALARS ! Scalars
  IF (LSCDERS) THEN
    KSCALARS_NSDER_OFFSET = IFIRST
    IFIRST = IFIRST + KF_SCALARS ! Scalars NS Derivatives
  ENDIF
  ! the rest of fields is being computed  in fourier space, namely in FSC
  IF_LEG = IFIRST
  IF (LUVDER) THEN
    KUV_EWDER_OFFSET = IFIRST
    IFIRST = IFIRST+2*KF_UV ! U and V derivatives
  ENDIF
  IF (LSCDERS) THEN
    KSCALARS_EWDER_OFFSET = IFIRST
    IFIRST = IFIRST + KF_SCALARS ! Scalars EW Derivatives
  ENDIF
  IF_FOURIER = IFIRST
  IF (IF_FOURIER /= KF_FS) CALL ABORT_TRANS('Size mismatch: Wrong computation KF_FS')

  ALLOCATOR = MAKE_BUFFERED_ALLOCATOR()
  HLTINV = PREPARE_LTINV(ALLOCATOR,KF_UV,KF_SCALARS,LVORGP,LDIVGP,LSCDERS)
  HTRMTOL = PREPARE_TRMTOL(ALLOCATOR,IF_LEG)
  HFOURIER_IN = PREPARE_FOURIER_IN(ALLOCATOR,IF_FOURIER)
  HFSC = PREPARE_FSC(ALLOCATOR)
  HFTINV = PREPARE_FTINV(ALLOCATOR)
  HTRLTOG = PREPARE_TRLTOG(ALLOCATOR,IF_FOURIER,KF_GP)

  CALL INSTANTIATE_ALLOCATOR(ALLOCATOR, REUSE_PTR)

  ! No splitting of fields, transform done in one go
  ! from PSPXXX to FOUBUF

  CALL GSTATS(102,0)
  IF (KF_FS > 0) THEN
    CALL LTINV(ALLOCATOR,HLTINV,KF_UV,KF_SCALARS,&
        & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2, &
        & FOUBUF_IN)
  ENDIF
  CALL GSTATS(102,1)

  CALL GSTATS(152,0)
  WRITE(NOUT,*) 'ltinv_ctl:TRMTOL_CUDAAWARE'
  IF (KF_FS > 0) CALL TRMTOL_CUDAAWARE(ALLOCATOR,HTRMTOL,FOUBUF_IN,FOUBUF,IF_LEG)
  CALL GSTATS(152,1)

  CALL GSTATS(107,0)
  ! from FOUBUF to PREEL_COMPLEX. Divide by two because we consider this complex space now
  IF (KF_FS > 0) CALL FOURIER_IN(ALLOCATOR,HFOURIER_IN,FOUBUF,PREEL_COMPLEX,IF_LEG,IF_FOURIER)

  !    2.  Fourier space computations

  ! fill the rest of PREEL_COMPLEX
  IF (KF_FS > 0) CALL FSC(ALLOCATOR,HFSC,PREEL_COMPLEX, IF_FOURIER, KF_UV, KF_SCALARS, KUV_OFFSET, KSCALARS_OFFSET, &
      & KSCALARS_NSDER_OFFSET, KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET)

  !   3.  Fourier transform
  ! inplace operation
  IF (KF_FS > 0) CALL FTINV(ALLOCATOR, HFTINV, PREEL_COMPLEX,PREEL_REAL,IF_FOURIER)

  CALL GSTATS(107,1)

  !   4.  Transposition

  IF (PRESENT(KVSETUV)) THEN
    IVSETUV(:) = KVSETUV(:)
  ELSE
    IVSETUV(:) = -1
  ENDIF
  IVSETSC(:)=-1
  IF (PRESENT(KVSETSC)) THEN
    IVSETSC(:) = KVSETSC(:)
  ELSE
    IOFF=0
    IF (PRESENT(KVSETSC2)) THEN
      IFGP2=UBOUND(KVSETSC2,1)
      IVSETSC(1:IFGP2)=KVSETSC2(:)
      IOFF=IOFF+IFGP2
    ENDIF
    IF (PRESENT(KVSETSC3A)) THEN
      IFGP3A=UBOUND(KVSETSC3A,1)
      IGP3APAR=UBOUND(PGP3A,3)
      IF (LSCDERS) IGP3APAR=IGP3APAR/3
      DO J3=1,IGP3APAR
        IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
        IOFF=IOFF+IFGP3A
      ENDDO
    ENDIF
    IF (PRESENT(KVSETSC3B)) THEN
      IFGP3B=UBOUND(KVSETSC3B,1)
      IGP3BPAR=UBOUND(PGP3B,3)
      IF (LSCDERS) IGP3BPAR=IGP3BPAR/3
      DO J3=1,IGP3BPAR
        IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
        IOFF=IOFF+IFGP3B
      ENDDO
    ENDIF
    IF (IOFF > 0 .AND. IOFF /= KF_SCALARS_G ) THEN
      WRITE(NERR,*)'FTINV:IOFF,KF_SCALARS_G ',IOFF,KF_SCALARS_G
      CALL ABORT_TRANS('FTINV_CTL_MOD:IOFF /= KF_SCALARS_G')
    ENDIF
  ENDIF

  IST = 1
  IF (KF_UV_G > 0) THEN
    IF (LVORGP) THEN
      IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
      IST = IST+KF_UV_G
    ENDIF
    IF ( LDIVGP) THEN
      IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
      IST = IST+KF_UV_G
    ENDIF
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IF (KF_SCALARS_G > 0) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
    IST = IST+KF_SCALARS_G
    IF (LSCDERS) THEN
      IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
      IST = IST+KF_SCALARS_G
    ENDIF
  ENDIF
  IF (KF_UV_G > 0 .AND. LUVDER) THEN
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
    IVSET(IST:IST+KF_UV_G-1) = IVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IF (KF_SCALARS_G > 0) THEN
    IF (LSCDERS) THEN
      IVSET(IST:IST+KF_SCALARS_G-1) = IVSETSC(:)
      IST = IST+KF_SCALARS_G
    ENDIF
  ENDIF

  CALL GSTATS(157,0)
  CALL TRLTOG_CUDAAWARE(ALLOCATOR,HTRLTOG,PREEL_REAL,IF_FOURIER,KF_GP,KF_UV_G,KF_SCALARS_G,IVSET,&
   &PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)
  CALL GSTATS(157,1)
 !     ------------------------------------------------------------------
!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_CTL
END MODULE INV_TRANS_CTL_MOD
