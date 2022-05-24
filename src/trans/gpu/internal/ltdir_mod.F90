! (C) Copyright 1987- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_MOD
  CONTAINS
  SUBROUTINE LTDIR(FOUBUF,KF_FS,KF_UV,KF_SCALARS,&
   & PSPVOR,PSPDIV,PSPSCALAR,&
   & PSPSC3A,PSPSC3B,PSPSC2, &
   & KFLDPTRUV,KFLDPTRSC)
  
  
  USE PARKIND_ECTRANS  ,ONLY : JPIM, JPRBT, JPRD, JPRB
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX
  
  USE TPM_DIM     ,ONLY : R
  USE TPM_DISTR   ,ONLY : D
  USE TPM_GEOMETRY
  
  USE PREPSNM_MOD ,ONLY : PREPSNM
  USE LEDIR_MOD   ,ONLY : LEDIR, LEDIR_POINTERS, LEDIR_PACK_BUFFER, LEDIR_ALLOC_SIZE
  USE UVTVD_MOD
  USE UPDSP_MOD   ,ONLY : UPDSP
  USE UPDSPB_MOD   ,ONLY : UPDSPB
  USE MPL_MODULE      ,ONLY : MPL_BARRIER
  USE TPM_GEN         ,ONLY : LSYNC_TRANS
  USE TPM_TRANS       ,ONLY : NF_SC2, NF_SC3A, NF_SC3B
  USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

USE TPM_TRANS, ONLY: REUSE_PTR
USE TPM_DISTR, ONLY : NPROC

  
  !**** *LTDIR* - Control of Direct Legendre transform step
  
  !     Purpose.
  !     --------
  !        Tranform from Fourier space to spectral space, compute
  !        vorticity and divergence.
  
  !**   Interface.
  !     ----------
  !        *CALL* *LTDIR(...)*
  
  !        Explicit arguments :
  !        --------------------  KM     - zonal wavenumber
  !                              KMLOC  - local zonal wavenumber
  
  !        Implicit arguments :  None
  !        --------------------
  
  !     Method.
  !     -------
  
  !     Externals.
  !     ----------
  !         PREPSNM - prepare REPSNM for wavenumber KM
  !         PRFI2   - prepares the Fourier work arrays for model variables.
  !         LEDIR   - direct Legendre transform
  !         UVTVD   -
  !         UPDSP   - updating of spectral arrays (fields)
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  
  !     Author.
  !     -------
  !        Mats Hamrud and Philippe Courtier  *ECMWF*
  
  !     Modifications.
  !     --------------
  !        Original : 87-11-24
  !        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
  !                            for uv formulation
  !        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
  !        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
  !        Modified 94-04-06 R. El khatib Full-POS implementation
  !        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
  !                             instead of u,v->vor,div
  !        MPP Group : 95-10-01 Support for Distributed Memory version
  !        K. YESSAD (AUGUST 1996):
  !               - Legendre transforms for transmission coefficients.
  !        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
  !        R. El Khatib 12-Jul-2012 LDSPC2 replaced by UVTVD
  !     ------------------------------------------------------------------
  
  IMPLICIT NONE
  
  !     DUMMY INTEGER SCALARS
  INTEGER(KIND=JPIM)  :: KM
  INTEGER(KIND=JPIM)  :: KMLOC
  INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS
  
  REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
  REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
  REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
  REAL(KIND=JPRBT),  ALLOCATABLE, INTENT(INOUT) :: FOUBUF(:)
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: IFC, IIFC, IDGLU, IFIRST
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  REAL(KIND=JPRB), ALLOCATABLE, TARGET :: POA1(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE, TARGET :: POA2(:,:,:)
  REAL(KIND=JPRB), POINTER :: PU(:,:,:), PV(:,:,:), PVOR(:,:,:), PDIV(:,:,:)

REAL(KIND=JPRBT), POINTER :: ZINPS(:), ZINPA(:)
REAL(KIND=JPRD), POINTER :: ZINPS0(:), ZINPA0(:)
REAL(KIND=JPRBT), POINTER :: ZOUT(:)
REAL(KIND=JPRD), POINTER :: ZOUT0(:)
INTEGER(KIND=8)  :: ALLOC_SZ, ALLOC_POS
  
  !     ------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',0,ZHOOK_HANDLE)
  
  !     ------------------------------------------------------------------
  
  !*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
  !              --------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  !*       2.    PREPARE WORK ARRAYS.
  !              --------------------
  IF (NPROC == 1) THEN
    ! Short cut - no need to go through tansforms, we will go directly into
    ! the legendre space, but for that we need twice the memory, roughly
    ! (but we don't need the send/recv buffers)
    ALLOC_SZ = KF_FS*D%NLENGTF*SIZEOF(REUSE_PTR(1))
  ELSE
    ALLOC_SZ = 0
  ENDIF
  
  ! Check if the reuse buffer is large enough
  ALLOC_SZ = ALLOC_SZ+LEDIR_ALLOC_SIZE(KF_FS)
  
  IF (.NOT. ALLOCATED(REUSE_PTR)) THEN
    ALLOCATE(REUSE_PTR(ALLOC_SZ/SIZEOF(REUSE_PTR(1))))
    !$ACC ENTER DATA CREATE(REUSE_PTR)
  ELSEIF (SIZEOF(REUSE_PTR) < ALLOC_SZ) THEN
    ! and reallocate if needed
    !$ACC EXIT DATA DELETE(REUSE_PTR)
    DEALLOCATE(REUSE_PTR)
    ALLOCATE(REUSE_PTR(ALLOC_SZ/SIZEOF(REUSE_PTR(1))))
    !$ACC ENTER DATA CREATE(REUSE_PTR)
  ENDIF
  
  ! Figure out which pointers to use
  ALLOC_POS=1
  IF (NPROC == 1) THEN
    ALLOC_POS = ALLOC_POS+KF_FS*D%NLENGTF
  ENDIF
  
  CALL LEDIR_POINTERS(REUSE_PTR(ALLOC_POS:),KF_FS,ZINPS,ZINPA,ZOUT,ZINPS0,ZINPA0,ZOUT0)
  
  CALL LEDIR_PACK_BUFFER(FOUBUF,ZINPS,ZINPA,ZINPS0,ZINPA0,KF_FS,KF_UV)
  IF (NPROC > 1) THEN
    !$ACC EXIT DATA DELETE(FOUBUF)
    DEALLOCATE(FOUBUF)
  ENDIF


  ALLOCATE(POA1(2*KF_FS,R%NTMAX+3,D%NUMP))
  !$ACC ENTER DATA CREATE(POA1)

  ! do the legendre transform
CALL LEDIR(ZINPS,ZINPA,ZOUT,ZINPS0,ZINPA0,ZOUT0,POA1,KF_FS)

  !$ACC DATA COPYOUT(PSPVOR,PSPDIV) IF(KF_UV > 0)
  !$ACC DATA COPYOUT(PSPSCALAR) IF(PRESENT(PSPSCALAR) .AND. KF_SCALARS > 0)
  !$ACC DATA COPYOUT(PSPSC2) IF(NF_SC2 > 0)
  !$ACC DATA COPYOUT(PSPSC3A) IF(NF_SC3A > 0)
  !$ACC DATA COPYOUT(PSPSC3B) IF(NF_SC3B > 0)

  !     ------------------------------------------------------------------
  
  !*       5.    COMPUTE VORTICITY AND DIVERGENCE.
  !              ---------------------------------

  IF( KF_UV > 0 ) THEN
     ALLOCATE(POA2(4*KF_UV,R%NTMAX+3,D%NUMP))
     !$ACC ENTER DATA CREATE(POA2)

     ! U and V are in POA1
     IFIRST = 0
     PU => POA1(IFIRST+1:IFIRST+2*KF_UV,:,:)
     IFIRST = IFIRST + 2*KF_UV
     PV => POA1(IFIRST+1:IFIRST+2*KF_UV,:,:)
     ! Compute VOR and DIV ino POA2
     IFIRST = 0
     PVOR => POA2(IFIRST+1:IFIRST+2*KF_UV,:,:)
     IFIRST = IFIRST + 2*KF_UV
     PDIV => POA2(IFIRST+1:IFIRST+2*KF_UV,:,:)

     ! Compute vorticity and divergence
     CALL UVTVD(KF_UV,PU,PV,PVOR,PDIV)

     ! Write back. Note, if we have UV, the contract says we *must* have VOR/DIV
     CALL UPDSPB(KF_UV,PVOR,PSPVOR,KFLDPTRUV)
     CALL UPDSPB(KF_UV,PDIV,PSPDIV,KFLDPTRUV)


     !$ACC EXIT DATA DELETE(POA2)
     DEALLOCATE(POA2)
  ENDIF
  !     ------------------------------------------------------------------
  
  !*       6.    UPDATE SPECTRAL ARRAYS.
  !              -----------------------
  
  ! this is on the host, so need to cp from device, Nils
  CALL UPDSP(KF_UV,KF_SCALARS,POA1,&
   & PSPSCALAR,&
   & PSPSC3A,PSPSC3B,PSPSC2 , &
   & KFLDPTRUV,KFLDPTRSC)
  
  !$ACC EXIT DATA DELETE(POA1)
  DEALLOCATE(POA1)

  IF (LSYNC_TRANS) THEN
    CALL GSTATS(430,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(430,1)
  ENDIF
  CALL GSTATS(412,0)
  !$ACC END DATA
  !$ACC END DATA
  !$ACC END DATA
  !$ACC END DATA
  !$ACC END DATA
  IF (LSYNC_TRANS) THEN
    CALL GSTATS(432,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(432,1)
  ENDIF
  CALL GSTATS(412,1)

  !     ------------------------------------------------------------------
  
  IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',1,ZHOOK_HANDLE)

  END SUBROUTINE LTDIR
  END MODULE LTDIR_MOD
