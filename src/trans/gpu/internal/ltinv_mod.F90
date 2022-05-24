! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV_MOD
  CONTAINS
  SUBROUTINE LTINV(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1,&
   & PSPVOR,PSPDIV,PSPSCALAR,&
   & PSPSC3A,PSPSC3B,PSPSC2 , &
   & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)
  
  USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,   JPRBT
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  
  USE TPM_DIM         ,ONLY : R
  USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B, FOUBUF_IN, LSCDERS
  USE TPM_FLT
  USE TPM_GEOMETRY
  USE TPM_DISTR       ,ONLY : D
  use tpm_gen, only: nout
  !USE PRLE1_MOD
  USE PREPSNM_MOD     ,ONLY : PREPSNM
  USE PRFI1B_MOD      ,ONLY : PRFI1B
  USE VDTUV_MOD       ,ONLY : VDTUV
  USE SPNSDE_MOD      ,ONLY : SPNSDE
  USE LEINV_MOD       ,ONLY : LEINV
  USE FSPGL_INT_MOD   ,ONLY : FSPGL_INT
  USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
  use ieee_arithmetic
  USE TPM_FIELDS      ,ONLY : F,ZEPSNM
  USE MPL_MODULE      ,ONLY : MPL_BARRIER
  USE TPM_GEN         ,ONLY : LSYNC_TRANS
  
  
  !**** *LTINV* - Inverse Legendre transform
  !
  !     Purpose.
  !     --------
  !        Tranform from Laplace space to Fourier space, compute U and V
  !        and north/south derivatives of state variables.
  
  !**   Interface.
  !     ----------
  !        *CALL* *LTINV(...)
  
  !        Explicit arguments :
  !        --------------------
  !          KM        - zonal wavenumber
  !          KMLOC     - local zonal wavenumber
  !          PSPVOR    - spectral vorticity
  !          PSPDIV    - spectral divergence
  !          PSPSCALAR - spectral scalar variables
  
  !        Implicit arguments :  The Laplace arrays of the model.
  !        --------------------  The values of the Legendre polynomials
  !                              The grid point arrays of the model
  !     Method.
  !     -------
  
  !     Externals.
  !     ----------
  
  !         PREPSNM - prepare REPSNM for wavenumber KM
  !         PRFI1B  - prepares the spectral fields
  !         VDTUV   - compute u and v from vorticity and divergence
  !         SPNSDE  - compute north-south derivatives
  !         LEINV   - Inverse Legendre transform
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  !        Temperton, 1991, MWR 119 p1303
  
  !     Author.
  !     -------
  !        Mats Hamrud  *ECMWF*
  
  !     Modifications.
  !     --------------
  !        Original : 00-02-01 From LTINV in IFS CY22R1
  !     ------------------------------------------------------------------
  
  IMPLICIT NONE
  
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2
  INTEGER(KIND=JPIM), INTENT(IN) :: KDIM1
  
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPVOR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPDIV(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSCALAR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC2(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3B(:,:,:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
  
  
  EXTERNAL  FSPGL_PROC
  OPTIONAL  FSPGL_PROC
  
  INTEGER(KIND=JPIM) :: IIFC, IDGLU
  INTEGER(KIND=JPIM) :: IFIRST, J3

  REAL(KIND=JPRB), ALLOCATABLE, TARGET :: PIA(:,:,:)
  REAL(KIND=JPRB), POINTER :: PU(:,:,:), PV(:,:,:), PVOR(:,:,:), PDIV(:,:,:)
  REAL(KIND=JPRB), POINTER :: PSCALARS(:,:,:), PSCALARS_NSDER(:,:,:)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  !CHARACTER(LEN=10) :: CLHOOK
  
  
  INTEGER(KIND=JPIM) :: KM
  INTEGER(KIND=JPIM) :: KMLOC
  
  
  !     ------------------------------------------------------------------
  
  !*       1.       PERFORM LEGENDRE TRANFORM.
  !                 --------------------------
  
  !WRITE(CLHOOK,FMT='(A,I4.4)') 'LTINV_',KM
  IF (LHOOK) CALL DR_HOOK('LTINV_MOD',0,ZHOOK_HANDLE)
  
  !     ------------------------------------------------------------------
  
  
  !*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
  !              ----------------------------------------------
  
! COPY FROM PSPXXXX TO ZIA

  IF (LSYNC_TRANS) THEN
    CALL MPL_BARRIER(CDSTRING='')
  ENDIF
  CALL GSTATS(431,0)
  !$ACC DATA COPYIN(PSPVOR,PSPDIV) IF(KF_UV > 0)
  !$ACC DATA COPYIN(PSPSCALAR) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSCALAR))
  !$ACC DATA COPYIN(PSPSC2) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC2) .AND. NF_SC2 > 0)
  !$ACC DATA COPYIN(PSPSC3A) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC3A) .AND. NF_SC3A > 0)
  !$ACC DATA COPYIN(PSPSC3B) IF(KF_SCALARS > 0 .AND. PRESENT(PSPSC3B) .AND. NF_SC3B > 0)
  IF (LSYNC_TRANS) THEN
    CALL MPL_BARRIER(CDSTRING='')
  ENDIF

  ! Compute PIA Domain decomposition
  IFIRST = 0
  IFIRST = IFIRST + 2*KF_UV ! Vorticity or divergence
  IFIRST = IFIRST + 2*KF_UV ! Vorticity or divergence
  IFIRST = IFIRST + 2*KF_UV ! U
  IFIRST = IFIRST + 2*KF_UV ! V
  IFIRST = IFIRST + 2*KF_SCALARS ! Scalars
  IF (LSCDERS) IFIRST = IFIRST + 2*KF_SCALARS ! Scalars NS Derivatives

  ! Allocate data accordingly
  ALLOCATE(PIA(IFIRST,R%NSMAX+3,D%NUMP))
  !$ACC ENTER DATA CREATE(PIA)

  ! And reiterate domain decomposition to assign pointers
  IFIRST = 0
  IF (.NOT. LVORGP .OR. LDIVGP) THEN
    ! Usually we want to store vorticity first
    PVOR => PIA(IFIRST+1:IFIRST+2*KF_UV,:,:)
    IFIRST = IFIRST + 2*KF_UV ! Vorticity

    PDIV => PIA(IFIRST+1:IFIRST+2*KF_UV,:,:)
    IFIRST = IFIRST + 2*KF_UV ! Divergence
  ELSE
    ! Except if we want to translate Vorticity but not Divergence, we should have Divergence first
    ! Then we have all buffers that move on in a contiguous buffer
    PDIV => PIA(IFIRST+1:IFIRST+2*KF_UV,:,:)
    IFIRST = IFIRST + 2*KF_UV ! Divergence

    PVOR => PIA(IFIRST+1:IFIRST+2*KF_UV,:,:)
    IFIRST = IFIRST + 2*KF_UV ! Vorticity
  ENDIF
  PU => PIA(IFIRST+1:IFIRST+2*KF_UV,:,:)
  IFIRST = IFIRST + 2*KF_UV ! U
  PV => PIA(IFIRST+1:IFIRST+2*KF_UV,:,:)
  IFIRST = IFIRST + 2*KF_UV ! V
  PSCALARS => PIA(IFIRST+1:IFIRST+2*KF_SCALARS,:,:)
  IFIRST = IFIRST + 2*KF_SCALARS ! Scalars
  PSCALARS_NSDER => PIA(IFIRST+1:IFIRST+2*KF_SCALARS,:,:)
  IF (LSCDERS) IFIRST = IFIRST + 2*KF_SCALARS ! Scalars NS Derivatives

  CALL GSTATS(431,1)
  IF (KF_UV > 0) THEN
    CALL PRFI1B(PVOR,PSPVOR,KF_UV,UBOUND(PSPVOR,2),KFLDPTRUV)
    CALL PRFI1B(PDIV,PSPDIV,KF_UV,UBOUND(PSPDIV,2),KFLDPTRUV)

    ! Compute U and V for VOR and DIV
    CALL VDTUV(KF_UV,ZEPSNM,PVOR,PDIV,PU,PV)
  ENDIF

  IF (KF_SCALARS > 0) THEN
    IF(PRESENT(PSPSCALAR)) THEN
      CALL PRFI1B(PSCALARS,PSPSCALAR,KF_SCALARS,UBOUND(PSPSCALAR,2),KFLDPTRSC)
    ELSE
      IFIRST = 1
      IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
        CALL PRFI1B(PSCALARS(IFIRST:IFIRST+2*NF_SC2-1,:,:),PSPSC2(:,:),NF_SC2,UBOUND(PSPSC2,2))
        IFIRST  = IFIRST+2*NF_SC2
      ENDIF
      IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
        DO J3=1,UBOUND(PSPSC3A,3)
          CALL PRFI1B(PSCALARS(IFIRST:IFIRST+2*NF_SC3A-1,:,:),PSPSC3A(:,:,J3),NF_SC3A,UBOUND(PSPSC3A,2))
          IFIRST  = IFIRST+2*NF_SC3A
        ENDDO
      ENDIF
      IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
        DO J3=1,UBOUND(PSPSC3B,3)
          CALL PRFI1B(PSCALARS(IFIRST:IFIRST+2*NF_SC3B-1,:,:),PSPSC3B(:,:,J3),NF_SC3B,UBOUND(PSPSC3B,2))
          IFIRST  = IFIRST+2*NF_SC3B
        ENDDO
      ENDIF
      IF(IFIRST-1 /= 2*KF_SCALARS) THEN
        WRITE(0,*) 'LTINV:KF_SCALARS,IFIRST',KF_SCALARS,IFIRST
        CALL ABORT_TRANS('LTINV_MOD:IFIRST /= 2*KF_SCALARS')
      ENDIF
    ENDIF
  ENDIF
  !$ACC END DATA
  !$ACC END DATA
  !$ACC END DATA
  !$ACC END DATA
  !$ACC END DATA

  IF (KF_SCDERS > 0) THEN
    CALL SPNSDE(KF_SCALARS,ZEPSNM,PSCALARS,PSCALARS_NSDER)
  ENDIF

  !     ------------------------------------------------------------------


  !*       4.    INVERSE LEGENDRE TRANSFORM.
  !              ---------------------------

  ! Forget aobut Vorticity and divergence if we don't need it in the output
  IFIRST = 1
  IF(.NOT. LVORGP) IFIRST = IFIRST+2*KF_UV
  IF(.NOT. LDIVGP) IFIRST = IFIRST+2*KF_UV

  ! Transform PIA into FOUBUF_IN
  IF (KF_OUT_LT > 0) THEN
    CALL LEINV(PIA(IFIRST:,:,:), FOUBUF_IN)
  ENDIF

  !$ACC EXIT DATA DELETE(PIA)
  DEALLOCATE(PIA)

  IF (LHOOK) CALL DR_HOOK('LTINV_MOD',1,ZHOOK_HANDLE)
  !     ------------------------------------------------------------------
  
  END SUBROUTINE LTINV
  END MODULE LTINV_MOD
  
