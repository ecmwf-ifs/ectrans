#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV_MOD
  USE BUFFERED_ALLOCATOR_MOD, ONLY: ALLOCATION_RESERVATION_HANDLE

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LTINV, LTINV_HANDLE, PREPARE_LTINV

  TYPE LTINV_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HPIA_AND_IN
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HOUTS_AND_OUTA
  END TYPE

CONTAINS
  FUNCTION PREPARE_LTINV(ALLOCATOR,KF_UV,KF_SCALARS,LVORGP,LDIVGP,LSCDERS) RESULT(HLTINV)
    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPRD
    USE TPM_DISTR,              ONLY: D
    USE TPM_DIM,                ONLY: R
    USE ISO_C_BINDING,          ONLY: C_SIZE_T, C_SIZEOF
    USE LEINV_MOD,              ONLY: LEINV_STRIDES
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, RESERVE

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV,KF_SCALARS
    LOGICAL, INTENT(IN) :: LVORGP,LDIVGP,LSCDERS

    TYPE(LTINV_HANDLE) :: HLTINV

    INTEGER(KIND=C_SIZE_T) :: IALLOC_SZ, IPIA_SZ
    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0, IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM)  :: IIN_STRIDES0, IIN_SIZE
    INTEGER(KIND=JPIM)  :: IIN0_STRIDES0, IIN0_SIZE

    REAL(KIND=JPRBT) :: ZPRBT_DUMMY
    REAL(KIND=JPRD) :: ZPRD_DUMMY

    INTEGER(KIND=JPIM) :: IF_READIN, IF_LEG

    ! # fields that are initially read. We always read vorticity
    ! and divergence! Also keep in mind that we actually have 2X
    ! this number of levels because real+complex
    IF_READIN = 0
    IF_READIN = IF_READIN + KF_UV ! Vorticity or divergence
    IF_READIN = IF_READIN + KF_UV ! Vorticity or divergence
    IF_READIN = IF_READIN + KF_UV ! U
    IF_READIN = IF_READIN + KF_UV ! V
    IF_READIN = IF_READIN + KF_SCALARS ! Scalars
    IF (LSCDERS) &
      IF_READIN = IF_READIN + KF_SCALARS ! Scalars NS Derivatives

    IPIA_SZ = ALIGN(INT(2*IF_READIN*(R%NSMAX+3)*D%NUMP,KIND=C_SIZE_T)*C_SIZEOF(ZPRBT_DUMMY),128)

    ! In Legendre space, we then ignore vorticity/divergence, if
    ! they don't need to be transformed.
    IF_LEG = IF_READIN
    IF(.NOT. LVORGP) IF_LEG = IF_LEG - KF_UV ! No vorticity needed
    IF(.NOT. LDIVGP) IF_LEG = IF_LEG - KF_UV ! No divergence needed

    CALL LEINV_STRIDES(IF_LEG,IOUT_STRIDES0,IOUT_SIZE,IIN_STRIDES0,IIN_SIZE,&
                       IOUT0_STRIDES0,IOUT0_SIZE,IIN0_STRIDES0,IIN0_SIZE)

    ! PIA
    IALLOC_SZ = IPIA_SZ
    ! ZINP
    IALLOC_SZ = IALLOC_SZ + ALIGN(INT(IIN_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZPRBT_DUMMY),128)
    ! ZINP0
    IALLOC_SZ = IALLOC_SZ + ALIGN(INT(IIN0_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZPRD_DUMMY),128)

    HLTINV%HPIA_AND_IN = RESERVE(ALLOCATOR, IALLOC_SZ)

    IALLOC_SZ = 0
    ! ZOUTA
    IALLOC_SZ = IALLOC_SZ + ALIGN(INT(IOUT_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZPRBT_DUMMY),128)
    ! ZOUTS
    IALLOC_SZ = IALLOC_SZ + ALIGN(INT(IOUT_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZPRBT_DUMMY),128)
    ! ZOUTA0
    IALLOC_SZ = IALLOC_SZ + ALIGN(INT(IOUT0_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZPRD_DUMMY),128)
    ! ZOUTS0
    IALLOC_SZ = IALLOC_SZ + ALIGN(INT(IOUT0_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZPRD_DUMMY),128)

    HLTINV%HOUTS_AND_OUTA = RESERVE(ALLOCATOR, IALLOC_SZ)

  END FUNCTION PREPARE_LTINV

  SUBROUTINE LTINV(ALLOCATOR,HLTINV,KF_UV,KF_SCALARS,&
     & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2, &
     & ZOUTS,ZOUTA,ZOUTS0,ZOUTA0)

    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRB, JPRBT, JPRD
    USE YOMHOOK,                ONLY: LHOOK, DR_HOOK, JPHOOK
    USE TPM_DIM,                ONLY: R
    USE TPM_TRANS,              ONLY: LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B, LSCDERS
    USE TPM_GEOMETRY,           ONLY: G
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, ASSIGN_PTR, GET_ALLOCATION
    USE TPM_DISTR,              ONLY: D
    USE PRFI1B_MOD,             ONLY: PRFI1B
    USE VDTUV_MOD,              ONLY: VDTUV
    USE SPNSDE_MOD,             ONLY: SPNSDE
    USE LEINV_MOD,              ONLY: LEINV_STRIDES, LEINV
    USE ABORT_TRANS_MOD,        ONLY: ABORT_TRANS
    USE TPM_FIELDS_GPU,         ONLY: FG
    USE MPL_MODULE,             ONLY: MPL_BARRIER,MPL_ALL_MS_COMM
    USE TPM_GEN,                ONLY: LSYNC_TRANS
    USE TPM_STATS,              ONLY: GSTATS => GSTATS_NVTX
    USE ISO_C_BINDING,          ONLY: C_SIZE_T, C_LOC, C_SIZEOF

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
    REAL(KIND=JPRBT), POINTER, INTENT(OUT) :: ZOUTS(:), ZOUTA(:)
    REAL(KIND=JPRD), POINTER, INTENT(OUT) :: ZOUTS0(:), ZOUTA0(:)

    INTEGER(KIND=JPIM) :: IFIRST, J3

    REAL(KIND=JPRB), POINTER :: PIA_L(:), PIA(:,:,:)
    REAL(KIND=JPRB), POINTER :: PU(:,:,:), PV(:,:,:), PVOR(:,:,:), PDIV(:,:,:)
    REAL(KIND=JPRB), POINTER :: PSCALARS(:,:,:), PSCALARS_NSDER(:,:,:)

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(LTINV_HANDLE), INTENT(IN) :: HLTINV

    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0, IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM)  :: IIN_STRIDES0, IIN_SIZE
    INTEGER(KIND=JPIM)  :: IIN0_STRIDES0, IIN0_SIZE

    INTEGER(KIND=JPIM) :: IF_READIN, IF_LEG
    INTEGER(KIND=C_SIZE_T) :: IALLOC_POS, IALLOC_SZ

    REAL(KIND=JPRBT), POINTER :: ZINP(:)
    REAL(KIND=JPRD), POINTER :: ZINP0(:)

    ASSOCIATE(ZEPSNM=>FG%ZEPSNM)

    !     ------------------------------------------------------------------

    !*       1.       PERFORM LEGENDRE TRANFORM.
    !                 --------------------------

    IF (LHOOK) CALL DR_HOOK('LTINV_MOD',0,ZHOOK_HANDLE)

    ! Get all pointers
    IF_READIN = 0
    IF_READIN = IF_READIN + KF_UV ! Vorticity or divergence
    IF_READIN = IF_READIN + KF_UV ! Vorticity or divergence
    IF_READIN = IF_READIN + KF_UV ! U
    IF_READIN = IF_READIN + KF_UV ! V
    IF_READIN = IF_READIN + KF_SCALARS ! Scalars
    IF (LSCDERS) &
      IF_READIN = IF_READIN + KF_SCALARS ! Scalars NS Derivatives

    ! In Legendre space, we then ignore vorticity/divergence, if
    ! they don't need to be transformed.
    IF_LEG = IF_READIN
    IF(.NOT. LVORGP) IF_LEG = IF_LEG - KF_UV ! No vorticity needed
    IF(.NOT. LDIVGP) IF_LEG = IF_LEG - KF_UV ! No divergence needed

    CALL LEINV_STRIDES(IF_LEG,IOUT_STRIDES0,IOUT_SIZE,IIN_STRIDES0,IIN_SIZE,&
                       IOUT0_STRIDES0,IOUT0_SIZE,IIN0_STRIDES0,IIN0_SIZE)

    IALLOC_POS = 1

    ! PIA
    IALLOC_SZ = ALIGN(INT(2*IF_READIN*(R%NTMAX+3)*D%NUMP,KIND=C_SIZE_T)*C_SIZEOF(PIA_L(1)),128)
    CALL ASSIGN_PTR(PIA_L, GET_ALLOCATION(ALLOCATOR, HLTINV%HPIA_AND_IN),&
        & IALLOC_POS, IALLOC_SZ)
    CALL C_F_POINTER(C_LOC(PIA_L), PIA, (/ 2*IF_READIN, R%NTMAX+3, D%NUMP /))
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZINP
    IALLOC_SZ = ALIGN(INT(IIN_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZINP(1)),128)
    CALL ASSIGN_PTR(ZINP, GET_ALLOCATION(ALLOCATOR, HLTINV%HPIA_AND_IN),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZINP0
    IALLOC_SZ = ALIGN(INT(IIN0_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZINP0(1)),128)
    CALL ASSIGN_PTR(ZINP0, GET_ALLOCATION(ALLOCATOR, HLTINV%HPIA_AND_IN),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    IALLOC_POS = 1

    ! ZOUTA
    IALLOC_SZ = ALIGN(INT(IOUT_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZOUTA(1)),128)
    CALL ASSIGN_PTR(ZOUTA, GET_ALLOCATION(ALLOCATOR, HLTINV%HOUTS_AND_OUTA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZOUTS
    IALLOC_SZ = ALIGN(INT(IOUT_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZOUTS(1)),128)
    CALL ASSIGN_PTR(ZOUTS, GET_ALLOCATION(ALLOCATOR, HLTINV%HOUTS_AND_OUTA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZOUTA0
    IALLOC_SZ = ALIGN(INT(IOUT0_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZOUTA0(1)),128)
    CALL ASSIGN_PTR(ZOUTA0, GET_ALLOCATION(ALLOCATOR, HLTINV%HOUTS_AND_OUTA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZOUTS0
    IALLOC_SZ = ALIGN(INT(IOUT0_SIZE,KIND=C_SIZE_T)*C_SIZEOF(ZOUTS0(1)),128)
    CALL ASSIGN_PTR(ZOUTS0, GET_ALLOCATION(ALLOCATOR, HLTINV%HOUTS_AND_OUTA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! Assign pointers do the different components of PIA
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

    !     ------------------------------------------------------------------


    !*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
    !              ----------------------------------------------

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(422,0)
#ifdef OMPGPU
    !$OMP TARGET DATA MAP(TO:PSPVOR,PSPDIV) IF(KF_UV > 0)
    !$OMP TARGET DATA MAP(TO:PSPSCALAR) IF(PRESENT(PSPSCALAR) .AND. KF_SCALARS > 0)
    !$OMP TARGET DATA MAP(TO:PSPSC2) IF(NF_SC2 > 0)
    !$OMP TARGET DATA MAP(TO:PSPSC3A) IF(NF_SC3A > 0)
    !$OMP TARGET DATA MAP(TO:PSPSC3B) IF(NF_SC3B > 0)
#endif
#ifdef ACCGPU
    !$ACC DATA COPYIN(PSPVOR,PSPDIV) IF(KF_UV > 0)
    !$ACC DATA COPYIN(PSPSCALAR) IF(PRESENT(PSPSCALAR) .AND. KF_SCALARS > 0)
    !$ACC DATA COPYIN(PSPSC2) IF(NF_SC2 > 0)
    !$ACC DATA COPYIN(PSPSC3A) IF(NF_SC3A > 0)
    !$ACC DATA COPYIN(PSPSC3B) IF(NF_SC3B > 0)
#endif
    IF (LSYNC_TRANS) THEN
      CALL GSTATS(442,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(442,1)
    ENDIF
    CALL GSTATS(422,1)

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

    ! Compute NS derivatives if needed
    IF (LSCDERS) THEN
      CALL SPNSDE(KF_SCALARS,ZEPSNM,PSCALARS,PSCALARS_NSDER)
    ENDIF

#ifdef OMPGPU
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
    !$ACC WAIT(1)
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
#endif

  !     ------------------------------------------------------------------


    !*       4.    INVERSE LEGENDRE TRANSFORM.
    !              ---------------------------

    ! Legendre transforms. When converting PIA into ZOUT, we ignore the first entries of LEINV.
    ! This is because vorticity and divergence is not necessarily converted to GP space.
    CALL LEINV(ALLOCATOR,PIA(2*(IF_READIN-IF_LEG)+1:IF_READIN,:,:),ZINP,ZINP0,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,IF_LEG)

    IF (LHOOK) CALL DR_HOOK('LTINV_MOD',1,ZHOOK_HANDLE)
    END ASSOCIATE
    !     ------------------------------------------------------------------
  END SUBROUTINE LTINV
END MODULE LTINV_MOD

