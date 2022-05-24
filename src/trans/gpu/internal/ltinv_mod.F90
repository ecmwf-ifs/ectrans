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
  USE ALLOCATOR_MOD

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LTINV, LTINV_HANDLE, PREPARE_LTINV

  TYPE LTINV_HANDLE
  END TYPE

  CONTAINS
  FUNCTION PREPARE_LTINV(ALLOCATOR, KF_FS, KF_UV) RESULT(HLTINV)
    USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT, JPRD
    USE TPM_DISTR, ONLY: D
    USE TPM_DIM, ONLY: R
    USE ISO_C_BINDING
    USE LEDIR_MOD
    USE ALLOCATOR_MOD

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS, KF_UV

    TYPE(LTINV_HANDLE) :: HLTINV
  END FUNCTION

  SUBROUTINE LTINV(ALLOCATOR,HLTINV,KF_UV,KF_SCALARS,&
   & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2, &
   & FOUBUF_IN,FOUBUF_KFIELD)

  USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,   JPRBT
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

  USE TPM_DIM         ,ONLY : R
  USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B, LSCDERS
  USE TPM_FLT
  USE TPM_GEOMETRY
  USE TPM_DISTR       ,ONLY : D
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
  USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX
  
  
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
  
  
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
  
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPVOR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPDIV(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSCALAR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC2(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3B(:,:,:)
  REAL(KIND=JPRB), ALLOCATABLE, INTENT(OUT)  :: FOUBUF_IN(:)
  INTEGER(KIND=JPIM), INTENT(OUT) :: FOUBUF_KFIELD

  INTEGER(KIND=JPIM) :: IFIRST, J3

  REAL(KIND=JPRB), ALLOCATABLE, TARGET :: PIA(:,:,:)
  REAL(KIND=JPRB), POINTER :: PU(:,:,:), PV(:,:,:), PVOR(:,:,:), PDIV(:,:,:)
  REAL(KIND=JPRB), POINTER :: PSCALARS(:,:,:), PSCALARS_NSDER(:,:,:)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
  TYPE(LTINV_HANDLE), INTENT(IN) :: HLTINV

  !     ------------------------------------------------------------------

  !*       1.       PERFORM LEGENDRE TRANFORM.
  !                 --------------------------

  IF (LHOOK) CALL DR_HOOK('LTINV_MOD',0,ZHOOK_HANDLE)

  !     ------------------------------------------------------------------


  !*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
  !              ----------------------------------------------

  IF (LSYNC_TRANS) THEN
    CALL GSTATS(440,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(440,1)
  ENDIF
  CALL GSTATS(422,0)
  !$ACC DATA COPYIN(PSPVOR,PSPDIV) IF(KF_UV > 0)
  !$ACC DATA COPYIN(PSPSCALAR) IF(PRESENT(PSPSCALAR) .AND. KF_SCALARS > 0)
  !$ACC DATA COPYIN(PSPSC2) IF(NF_SC2 > 0)
  !$ACC DATA COPYIN(PSPSC3A) IF(NF_SC3A > 0)
  !$ACC DATA COPYIN(PSPSC3B) IF(NF_SC3B > 0)
  IF (LSYNC_TRANS) THEN
    CALL GSTATS(442,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(442,1)
  ENDIF
  CALL GSTATS(422,1)

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
  IF (LSCDERS) THEN
    PSCALARS_NSDER => PIA(IFIRST+1:IFIRST+2*KF_SCALARS,:,:)
    IFIRST = IFIRST + 2*KF_SCALARS ! Scalars NS Derivatives
  ENDIF

  IF (KF_UV > 0) THEN
    CALL PRFI1B(PVOR,PSPVOR,KF_UV,UBOUND(PSPVOR,2))
    CALL PRFI1B(PDIV,PSPDIV,KF_UV,UBOUND(PSPDIV,2))

    ! Compute U and V for VOR and DIV
    CALL VDTUV(KF_UV,ZEPSNM,PVOR,PDIV,PU,PV)
  ENDIF

  IF (KF_SCALARS > 0) THEN
    IF(PRESENT(PSPSCALAR)) THEN
      CALL PRFI1B(PSCALARS,PSPSCALAR,KF_SCALARS,UBOUND(PSPSCALAR,2))
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

  IF (LSCDERS) THEN
    CALL SPNSDE(KF_SCALARS,ZEPSNM,PSCALARS,PSCALARS_NSDER)
  ENDIF

  !     ------------------------------------------------------------------


  !*       4.    INVERSE LEGENDRE TRANSFORM.
  !              ---------------------------

  ! Forget about Vorticity and divergence if we don't need it in the output
  IFIRST = 1
  IF(.NOT. LVORGP) IFIRST = IFIRST+2*KF_UV
  IF(.NOT. LDIVGP) IFIRST = IFIRST+2*KF_UV

  ! Keep this for next functions because we have to remember this
  FOUBUF_KFIELD = SIZE(PIA,1)-IFIRST+1

  ! Transform PIA into FOUBUF_IN
  IF (FOUBUF_KFIELD > 0) THEN
    CALL LEINV(PIA(IFIRST:,:,:), FOUBUF_IN)
  ENDIF

  !$ACC EXIT DATA DELETE(PIA)
  DEALLOCATE(PIA)

  IF (LHOOK) CALL DR_HOOK('LTINV_MOD',1,ZHOOK_HANDLE)
  !     ------------------------------------------------------------------
  
  END SUBROUTINE LTINV
  END MODULE LTINV_MOD
  
