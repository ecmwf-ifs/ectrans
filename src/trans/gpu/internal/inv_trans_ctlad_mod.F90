! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS_CTLAD_MOD
  CONTAINS
    SUBROUTINE INV_TRANS_CTLAD(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT,&
      & KF_UV,KF_SCALARS,KF_SCDERS,&
      & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
      & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)
  
      !**** *INV_TRANS_CTLAD* - Control routine for adjoint of the inverse spectral transform.
  
      !     Purpose.
      !     --------
      !        Control routine for the adjoint of the inverse spectral transform
  
      !**   Interface.
      !     ----------
      !     CALL INV_TRANS_CTLAD(...)
  
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
  
      !     Externals.  LTINV_CTL   - control of Legendre transform
      !     ----------  FTINV_CTL   - control of Fourier transform
  
      !     Author.
      !     -------
      !        Mats Hamrud *ECMWF*
  
      !     Modifications.
      !     --------------
      !        Original : 01-01-03
  
      !     ------------------------------------------------------------------
  
  
      USE PARKIND_ECTRANS,        ONLY: JPIM, JPRB, JPRBT, JPRD
      USE ISO_C_BINDING,          ONLY: C_INT8_T
      USE TPM_GEN,                ONLY: NOUT
      USE TPM_TRANS,              ONLY: LDIVGP, LSCDERS, LUVDER, LVORGP, GROWING_ALLOCATION
      USE ABORT_TRANS_MOD,        ONLY: ABORT_TRANS
      USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, MAKE_BUFFERED_ALLOCATOR, INSTANTIATE_ALLOCATOR
      USE TRMTOLAD_MOD,           ONLY: PREPARE_TRMTOLAD, TRMTOLAD_HANDLE, TRMTOLAD
      USE LTINVAD_MOD,            ONLY: PREPARE_LTINVAD, LTINVAD_HANDLE, LTINVAD
      USE TRMTOLAD_PACK_UNPACK,   ONLY: TRMTOLAD_PACK_HANDLE, TRMTOLAD_UNPACK_HANDLE, &
        &                               PREPARE_TRMTOLAD_PACK, PREPARE_TRMTOLAD_UNPACK, TRMTOLAD_PACK, &
        &                               TRMTOLAD_UNPACK
      USE FSCAD_MOD,              ONLY: FSCAD_HANDLE, PREPARE_FSCAD, FSCAD
      USE FTDIR_MOD,              ONLY: FTDIR_HANDLE, PREPARE_FTDIR, FTDIR
      USE TRGTOL_MOD,             ONLY: TRGTOL_HANDLE, PREPARE_TRGTOL, TRGTOL
  
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
      REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
      REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
      REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
      REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
      REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
      REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
      INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
      INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
      INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
      INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
      INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
      REAL(KIND=JPRB) ,OPTIONAL   ,INTENT(IN) :: PGP(:,:,:)
      REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGPUV(:,:,:,:)
      REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP3A(:,:,:,:)
      REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP3B(:,:,:,:)
      REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP2(:,:,:)
  
      ! Local variables
  
      REAL(KIND=JPRB), POINTER :: FOUBUF(:), FOUBUF_IN(:)
      REAL(KIND=JPRBT), POINTER :: PREEL_REAL(:), PREEL_COMPLEX(:)
      REAL(KIND=JPRBT), POINTER :: ZOUTS(:), ZOUTA(:)
      REAL(KIND=JPRD), POINTER :: ZOUTS0(:), ZOUTA0(:)
      INTEGER(KIND=JPIM) :: KUV_OFFSET, KSCALARS_OFFSET, KSCALARS_NSDER_OFFSET, &
          & KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET
      INTEGER(KIND=JPIM) :: IF_LEG, IF_FOURIER
  
      INTEGER(KIND=JPIM) :: IFIRST
  
      TYPE(BUFFERED_ALLOCATOR) :: ALLOCATOR
      TYPE(LTINVAD_HANDLE) :: HLTINVAD
      TYPE(TRMTOLAD_PACK_HANDLE) :: HTRMTOL_PACK
      TYPE(TRMTOLAD_HANDLE) :: HTRMTOL
      TYPE(TRMTOLAD_UNPACK_HANDLE) :: HTRMTOL_UNPACK
      TYPE(FSCAD_HANDLE) :: HFSC
      TYPE(FTDIR_HANDLE) :: HFTDIR
      TYPE(TRGTOL_HANDLE) :: HTRGTOL
  
      INTEGER(KIND=C_INT8_T), POINTER :: PTR(:)
  
      !     ------------------------------------------------------------------

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
      HTRGTOL = PREPARE_TRGTOL(ALLOCATOR,KF_GP,IF_FOURIER)
      IF (KF_FS > 0) THEN
        HFTDIR = PREPARE_FTDIR(ALLOCATOR,IF_FOURIER)
        HFSC = PREPARE_FSCAD(ALLOCATOR)
        HTRMTOL_UNPACK = PREPARE_TRMTOLAD_UNPACK(ALLOCATOR,IF_LEG)
        HTRMTOL = PREPARE_TRMTOLAD(ALLOCATOR,IF_LEG)
        HTRMTOL_PACK = PREPARE_TRMTOLAD_PACK(ALLOCATOR,KF_UV,KF_SCALARS,LVORGP,LDIVGP,LSCDERS)
        HLTINVAD = PREPARE_LTINVAD(ALLOCATOR,KF_UV,KF_SCALARS,LVORGP,LDIVGP,LSCDERS)
      ENDIF
  
      CALL INSTANTIATE_ALLOCATOR(ALLOCATOR, GROWING_ALLOCATION)
  
      ! Adjoint of transposition into grid-point space
      CALL GSTATS(157,0)
      CALL TRGTOL(ALLOCATOR,HTRGTOL,PREEL_REAL,IF_FOURIER,KF_GP,KF_UV_G,KF_SCALARS_G,&
       & KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
       & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
       & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)
      CALL GSTATS(157,1)
  
      IF (KF_FS > 0) THEN
        CALL GSTATS(107,0)
        ! Fourier transformations
        CALL FTDIR(ALLOCATOR,HFTDIR,PREEL_REAL,PREEL_COMPLEX,IF_FOURIER)
        ! compute NS derivatives
        CALL FSCAD(ALLOCATOR,HFSC,PREEL_COMPLEX, IF_FOURIER, KF_UV, KF_SCALARS, KUV_OFFSET, KSCALARS_OFFSET, &
        & KSCALARS_NSDER_OFFSET, KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET)
        CALL GSTATS(107,1)
        
        ! Packing into send buffer, to fourier space and unpack
        CALL GSTATS(152,0)
        CALL TRMTOLAD_UNPACK(ALLOCATOR,HTRMTOL_UNPACK,FOUBUF,PREEL_COMPLEX,IF_LEG,IF_FOURIER)
        CALL TRMTOLAD(ALLOCATOR,HTRMTOL,FOUBUF_IN,FOUBUF,IF_LEG)
        CALL TRMTOLAD_PACK(ALLOCATOR,HTRMTOL_PACK,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,FOUBUF_IN,IF_LEG)
        CALL GSTATS(152,1)
        
        ! Legendre transformations
        CALL GSTATS(102,0)
        CALL LTINVAD(ALLOCATOR,HLTINVAD,KF_UV,KF_SCALARS,&
            & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2, &
            & ZOUTS,ZOUTA,ZOUTS0,ZOUTA0)
        CALL GSTATS(102,1)
      ENDIF
  
    END SUBROUTINE INV_TRANS_CTLAD
  END MODULE INV_TRANS_CTLAD_MOD
  