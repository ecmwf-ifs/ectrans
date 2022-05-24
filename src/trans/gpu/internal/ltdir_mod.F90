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
  SUBROUTINE LTDIR(KF_FS,KF_UV,KF_SCALARS,&
   & PSPVOR,PSPDIV,PSPSCALAR,&
   & PSPSC3A,PSPSC3B,PSPSC2, &
   & KFLDPTRUV,KFLDPTRSC)
  
  
  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  
  USE TPM_DIM     ,ONLY : R
  USE TPM_DISTR   ,ONLY : D
  USE TPM_GEOMETRY
  
  USE PREPSNM_MOD ,ONLY : PREPSNM
  USE LEDIR_MOD   ,ONLY : LEDIR
  USE UVTVD_MOD
  USE UPDSP_MOD   ,ONLY : UPDSP
  USE UPDSPB_MOD   ,ONLY : UPDSPB
   
  USE TPM_FIELDS      ,ONLY : ZOA1,ZOA2,ZEPSNM
  
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
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: IFC, IIFC, IDGLU
  INTEGER(KIND=JPIM) :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE
  
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  REAL(KIND=JPRB), ALLOCATABLE :: POA2(:,:,:)
  
  
  !call cudaProfilerStart
  
  !     ------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',0,ZHOOK_HANDLE)
  
  !     ------------------------------------------------------------------
  
  !*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
  !              --------------------------------------
  
  
  !     ------------------------------------------------------------------
  
  !*       2.    PREPARE WORK ARRAYS.
  !              --------------------
  
  ! do the legendre transform
  CALL LEDIR(KF_FS,KF_UV,ZOA1)

  !     ------------------------------------------------------------------
  
  !*       5.    COMPUTE VORTICITY AND DIVERGENCE.
  !              ---------------------------------

  IF( KF_UV > 0 ) THEN
     !!CALL PREPSNM

     IUS = 1
     IUE = 2*KF_UV
     IVS = 2*KF_UV+1
     IVE = 4*KF_UV
     IVORS = 1
     IVORE = 2*KF_UV
     IDIVS = 2*KF_UV+1
     IDIVE = 4*KF_UV

     ALLOCATE(POA2(4*KF_UV,R%NTMAX+3,D%NUMP))
     !$ACC ENTER DATA CREATE(POA2)


     ! Compute vorticity and divergence
     CALL UVTVD(KF_UV,ZOA1(IUS:IUE,:,:),ZOA1(IVS:IVE,:,:),&
         & POA2(IVORS:IVORE,:,:),POA2(IDIVS:IDIVE,:,:))

     ! Write back. Note, if we have UV, the contract says we *must* have VOR/DIV
     CALL UPDSPB(KF_UV,POA2(IVORS:IVORE,:,:),PSPVOR,KFLDPTRUV)
     CALL UPDSPB(KF_UV,POA2(IDIVS:IDIVE,:,:),PSPDIV,KFLDPTRUV)


     !$ACC EXIT DATA DELETE(POA2)
     DEALLOCATE(POA2)
  ENDIF
  !     ------------------------------------------------------------------
  
  !*       6.    UPDATE SPECTRAL ARRAYS.
  !              -----------------------
  
  !end loop over wavenumber
  
  !END DO
  
  !loop over wavenumber
  !DO KMLOC=1,D%NUMP
  !    KM = D%MYMS(KMLOC)
 
  ! this is on the host, so need to cp from device, Nils
  CALL UPDSP(KF_UV,KF_SCALARS,ZOA1,&
   & PSPSCALAR,&
   & PSPSC3A,PSPSC3B,PSPSC2 , &
   & KFLDPTRUV,KFLDPTRSC)
  
  !     ------------------------------------------------------------------
  
  IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',1,ZHOOK_HANDLE)
  
  !end loop over wavenumber
  !END DO

  
  !call cudaProfilerStop
  END SUBROUTINE LTDIR
  END MODULE LTDIR_MOD
