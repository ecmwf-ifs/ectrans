#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
! (C) Copyright 1987- ECMWF.
! (C) Copyright 1987- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIRAD_MOD
  USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPRB, JPRD, JPIB
  USE BUFFERED_ALLOCATOR_MOD, ONLY: ALLOCATION_RESERVATION_HANDLE
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PREPARE_LTDIRAD, LTDIRAD_HANDLE, LTDIRAD

  TYPE LTDIRAD_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HOUT_AND_POA
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HINPS_AND_ZINPA
  END TYPE

CONTAINS
  FUNCTION PREPARE_LTDIRAD(ALLOCATOR, KF_FS, KF_UV) RESULT(HLTDIR)
    USE TPM_DISTR,              ONLY: D
    USE TPM_DIM,                ONLY: R
    USE ISO_C_BINDING,          ONLY: C_SIZE_T, C_SIZEOF
    USE LEDIR_MOD,              ONLY: LEDIR_STRIDES
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, RESERVE

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS, KF_UV
    TYPE(LTDIRAD_HANDLE) :: HLTDIR

    INTEGER(KIND=JPIB) :: IALLOC_SZ
    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0
    INTEGER(KIND=JPIB)  :: IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM)  :: IIN_STRIDES0
    INTEGER(KIND=JPIB)  :: IIN_SIZE
    INTEGER(KIND=JPIM)  :: IIN0_STRIDES0, IIN0_SIZE

    REAL(KIND=JPRBT) :: ZPRBT_DUMMY
    REAL(KIND=JPRD) :: ZPRD_DUMMY

    CALL LEDIR_STRIDES(KF_FS,IOUT_STRIDES0=IOUT_STRIDES0,IOUT_SIZE=IOUT_SIZE,&
                       IOUT0_STRIDES0=IOUT0_STRIDES0,IOUT0_SIZE=IOUT0_SIZE,&
                       IIN_STRIDES0=IIN_STRIDES0,IIN_SIZE=IIN_SIZE,&
                       IIN0_STRIDES0=IIN0_STRIDES0,IIN0_SIZE=IIN0_SIZE)

    ! POA1
    IALLOC_SZ = ALIGN(2_JPIB*KF_FS*(R%NTMAX+3)*D%NUMP*C_SIZEOF(ZPRBT_DUMMY),128)
    ! POA2
    IALLOC_SZ = IALLOC_SZ + ALIGN(4_JPIB*KF_UV*(R%NTMAX+3)*D%NUMP*C_SIZEOF(ZPRBT_DUMMY),128)
    ! ZOUT
    IALLOC_SZ = IALLOC_SZ + ALIGN(IOUT_SIZE*C_SIZEOF(ZPRBT_DUMMY),128)
    ! ZOUT0
    IALLOC_SZ = IALLOC_SZ+ ALIGN(IOUT0_SIZE*C_SIZEOF(ZPRD_DUMMY),128)

    HLTDIR%HOUT_AND_POA = RESERVE(ALLOCATOR, IALLOC_SZ, "HLTDIRAD%HOUT_AND_POA")

    ! Check if the reuse buffer is large enough
    IALLOC_SZ = ALIGN(IIN_SIZE*C_SIZEOF(ZPRBT_DUMMY),128)
    IALLOC_SZ = IALLOC_SZ + ALIGN(IIN_SIZE*C_SIZEOF(ZPRBT_DUMMY),128)
    IALLOC_SZ = IALLOC_SZ + ALIGN(IIN0_SIZE*C_SIZEOF(ZPRD_DUMMY),128)
    IALLOC_SZ = IALLOC_SZ + ALIGN(IIN0_SIZE*C_SIZEOF(ZPRD_DUMMY),128)

    HLTDIR%HINPS_AND_ZINPA = RESERVE(ALLOCATOR, IALLOC_SZ, "HLTDIRAD%HINPS_AND_ZINPA")
  END FUNCTION PREPARE_LTDIRAD

  SUBROUTINE LTDIRAD(ALLOCATOR,HLTDIR,ZINPS,ZINPA,ZINPS0,ZINPA0,KF_FS,KF_UV,KF_SCALARS,&
                   & PSPVOR,PSPDIV,PSPSCALAR,&
                   & PSPSC3A,PSPSC3B,PSPSC2, &
                   & KFLDPTRUV,KFLDPTRSC)

    USE YOMHOOK,                ONLY: LHOOK, DR_HOOK, JPHOOK
    USE TPM_DIM,                ONLY: R
    USE TPM_DISTR,              ONLY: D
    USE TPM_GEOMETRY,           ONLY: G
    USE PREPSNM_MOD,            ONLY: PREPSNM
    USE LEDIR_MOD,              ONLY: LEDIR_STRIDES
    USE LEINV_MOD,              ONLY: LEINV
    USE UVTVDAD_MOD,            ONLY: UVTVDAD
    USE UPDSPAD_MOD,            ONLY: UPDSPAD
    USE UPDSPBAD_MOD,           ONLY: UPDSPBAD
    USE MPL_MODULE,             ONLY: MPL_BARRIER,MPL_ALL_MS_COMM
    USE TPM_GEN,                ONLY: LSYNC_TRANS, NERR
    USE TPM_TRANS,              ONLY: NF_SC2, NF_SC3A, NF_SC3B
    USE TPM_STATS,              ONLY: GSTATS => GSTATS_NVTX
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, ASSIGN_PTR, GET_ALLOCATION
    USE ISO_C_BINDING,          ONLY: C_SIZE_T, C_F_POINTER, C_LOC, C_SIZEOF

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

    REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPVOR(:,:)
    REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPDIV(:,:)
    REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSCALAR(:,:)
    REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC2(:,:)
    REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3A(:,:,:)
    REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3B(:,:,:)
    INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
    INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
    REAL(KIND=JPRBT)  ,POINTER, INTENT(OUT) :: ZINPS(:), ZINPA(:)
    REAL(KIND=JPRD)   ,POINTER, INTENT(OUT) :: ZINPS0(:), ZINPA0(:)

    !     LOCAL INTEGER SCALARS
    INTEGER(KIND=JPIM) :: IFC, IIFC, IDGLU, IFIRST

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JPRB), POINTER :: POA1_L(:), POA1(:,:,:)
    REAL(KIND=JPRB), POINTER :: POA2_L(:), POA2(:,:,:)
    REAL(KIND=JPRB), POINTER :: PU(:,:,:), PV(:,:,:), PVOR(:,:,:), PDIV(:,:,:)
    REAL(KIND=JPRBT), POINTER :: ZOUT(:)
    REAL(KIND=JPRD), POINTER :: ZOUT0(:)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(LTDIRAD_HANDLE), INTENT(IN) :: HLTDIR
    INTEGER(KIND=JPIB) :: IALLOC_POS, IALLOC_SZ
    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0
    INTEGER(KIND=JPIB)  :: IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM)  :: IIN_STRIDES0
    INTEGER(KIND=JPIB)  :: IIN_SIZE
    INTEGER(KIND=JPIM)  :: IIN0_STRIDES0, IIN0_SIZE

    INTEGER(KIND=JPIM)  :: JFLD, JN
  
    ASSOCIATE(D_NUMP=>D%NUMP, R_NTMAX=>R%NTMAX)

    !     ------------------------------------------------------------------
    IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',0,ZHOOK_HANDLE)

    !     ------------------------------------------------------------------

    !*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
    !              --------------------------------------


    !     ------------------------------------------------------------------

    !*       2.    PREPARE WORK ARRAYS.
    !              --------------------
    CALL LEDIR_STRIDES(KF_FS,IOUT_STRIDES0=IOUT_STRIDES0,IOUT_SIZE=IOUT_SIZE,&
                       IOUT0_STRIDES0=IOUT0_STRIDES0,IOUT0_SIZE=IOUT0_SIZE,&
                       IIN_STRIDES0=IIN_STRIDES0,IIN_SIZE=IIN_SIZE,&
                       IIN0_STRIDES0=IIN0_STRIDES0,IIN0_SIZE=IIN0_SIZE)

    IALLOC_POS = 1

    IALLOC_SZ = ALIGN(2_JPIB*KF_FS*(R%NTMAX+3)*D%NUMP*C_SIZEOF(POA1_L(1)),128)
    CALL ASSIGN_PTR(POA1_L, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ, SET_STREAM=1)
    CALL C_F_POINTER(C_LOC(POA1_L), POA1, (/ 2*KF_FS, R%NTMAX+3, D%NUMP /))
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    IALLOC_SZ = ALIGN(4_JPIB*KF_UV*(R%NTMAX+3)*D%NUMP*C_SIZEOF(POA2_L(1)),128)
    CALL ASSIGN_PTR(POA2_L, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ, SET_STREAM=1)
    CALL C_F_POINTER(C_LOC(POA2_L), POA2, (/ 4*KF_UV, R%NTMAX+3, D%NUMP /))
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZOUT
    IALLOC_SZ = ALIGN(IOUT_SIZE*C_SIZEOF(ZOUT(1)),128)
    CALL ASSIGN_PTR(ZOUT, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ, SET_STREAM=1)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZOUT0
    IALLOC_SZ = ALIGN(IOUT0_SIZE*C_SIZEOF(ZOUT0(1)),128)
    CALL ASSIGN_PTR(ZOUT0, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ, SET_STREAM=1)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    IALLOC_POS=1

    IALLOC_SZ = ALIGN(IIN_SIZE*C_SIZEOF(ZINPS(0)),128)
    CALL ASSIGN_PTR(ZINPS, GET_ALLOCATION(ALLOCATOR, HLTDIR%HINPS_AND_ZINPA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS=IALLOC_POS+IALLOC_SZ

    IALLOC_SZ = ALIGN(IIN_SIZE*C_SIZEOF(ZINPA(0)),128)
    CALL ASSIGN_PTR(ZINPA, GET_ALLOCATION(ALLOCATOR, HLTDIR%HINPS_AND_ZINPA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS=IALLOC_POS+IALLOC_SZ

    IALLOC_SZ = ALIGN(IIN0_SIZE*C_SIZEOF(ZINPS0(0)),128)
    CALL ASSIGN_PTR(ZINPS0, GET_ALLOCATION(ALLOCATOR, HLTDIR%HINPS_AND_ZINPA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS=IALLOC_POS+IALLOC_SZ

    IALLOC_SZ = ALIGN(IIN0_SIZE*C_SIZEOF(ZINPA0(0)),128)
    CALL ASSIGN_PTR(ZINPA0, GET_ALLOCATION(ALLOCATOR, HLTDIR%HINPS_AND_ZINPA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS=IALLOC_POS+IALLOC_SZ

#ifdef ACCGPU
    !$ACC DATA PRESENT(POA1,R_NTMAX,D_NUMP) 
    !$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(NONE) COPYIN(KF_FS)
#endif
#ifdef OMPGPU
    !$OMP TARGET DATA MAP(PRESENT,ALLOC:POA1,R_NTMAX,D_NUMP)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) DEFAULT(NONE) MAP(TO:KF_FS)
#endif
    DO KMLOC=1,D_NUMP
        DO JN=1,R_NTMAX+3
          DO JFLD=1,2*KF_FS
            POA1(JFLD,JN,KMLOC) = 0
          END DO
        END DO
    END DO
#ifdef ACCGPU
    !$ACC END DATA
#endif
#ifdef OMPGPU
    !$OMP END TARGET DATA
#endif

    !     ------------------------------------------------------------------

    !*       3.    PREPARE FOURIER ARRAYS.
    !              ----------------------

    !     ------------------------------------------------------------------

    !*       4.    COPY WORK ARRAYS TO DEVICE.
    !              ---------------------------


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
    !     ------------------------------------------------------------------

    !*       6.    UPDATE SPECTRAL ARRAYS.
    !              -----------------------

    ! this is on the host, so need to cp from device, Nils
    CALL UPDSPAD(KF_UV,KF_SCALARS,POA1,&
               & PSPSCALAR,&
               & PSPSC3A,PSPSC3B,PSPSC2 , &
               & KFLDPTRUV,KFLDPTRSC)

    !     ------------------------------------------------------------------

    !*       5.    COMPUTE VORTICITY AND DIVERGENCE.
    !              ---------------------------------

    IF( KF_UV > 0 ) THEN
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

       ! Write back. Note, if we have UV, the contract says we *must* have VOR/DIV
       CALL UPDSPBAD(KF_UV,PVOR,PSPVOR,KFLDPTRUV)
       CALL UPDSPBAD(KF_UV,PDIV,PSPDIV,KFLDPTRUV)
       
       ! Compute vorticity and divergence
       CALL UVTVDAD(KF_UV,PU,PV,PVOR,PDIV)
    ENDIF


#ifdef ACCGPU
    !$ACC WAIT(1)
#endif

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(430,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(430,1)
    ENDIF
    CALL GSTATS(412,0)
#ifdef OMPGPU
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
    !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
#endif

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(432,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(432,1)
    ENDIF
    CALL GSTATS(412,1)

    ! do the legendre transform
    CALL LEINV(ALLOCATOR,POA1,ZOUT,ZOUT0,ZINPS,ZINPA,ZINPS0,ZINPA0,KF_FS)

    !     ------------------------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('LTDIRAD_MOD',1,ZHOOK_HANDLE)

    END ASSOCIATE

  END SUBROUTINE LTDIRAD
END MODULE LTDIRAD_MOD
