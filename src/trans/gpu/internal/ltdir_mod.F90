#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
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
  USE ALLOCATOR_MOD
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PREPARE_LTDIR, LTDIR_HANDLE, LTDIR

  TYPE LTDIR_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HOUT_AND_POA
  END TYPE

CONTAINS
  FUNCTION PREPARE_LTDIR(ALLOCATOR, KF_FS, KF_UV) RESULT(HLTDIR)
    USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT, JPRD
    USE TPM_DISTR, ONLY: D
    USE TPM_DIM, ONLY: R
    USE ISO_C_BINDING
    USE LEDIR_MOD
    USE ALLOCATOR_MOD

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS, KF_UV
    TYPE(LTDIR_HANDLE) :: HLTDIR

    INTEGER(KIND=C_SIZE_T) :: IALLOC_SZ
    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0, IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE

    REAL(KIND=JPRBT) :: ZPRBT_DUMMY
    REAL(KIND=JPRD) :: ZPRD_DUMMY

    CALL LEDIR_STRIDES(KF_FS,IOUT_STRIDES0=IOUT_STRIDES0,IOUT_SIZE=IOUT_SIZE,&
                     IOUT0_STRIDES0=IOUT0_STRIDES0,IOUT0_SIZE=IOUT0_SIZE)

    ! POA1
    IALLOC_SZ = ALIGN(2*KF_FS*(R%NTMAX+3)*D%NUMP*SIZEOF(ZPRBT_DUMMY),128)
    ! POA2
    IALLOC_SZ = IALLOC_SZ + ALIGN(4*KF_UV*(R%NTMAX+3)*D%NUMP*SIZEOF(ZPRBT_DUMMY),128)
    ! ZOUT
    IALLOC_SZ = IALLOC_SZ + ALIGN(IOUT_SIZE*SIZEOF(ZPRBT_DUMMY),128)
    ! ZOUT0
    IALLOC_SZ = IALLOC_SZ+ ALIGN(IOUT0_SIZE*SIZEOF(ZPRD_DUMMY),128)

    HLTDIR%HOUT_AND_POA = RESERVE(ALLOCATOR, IALLOC_SZ)
  END FUNCTION

  SUBROUTINE LTDIR(ALLOCATOR,HLTDIR,ZINPS,ZINPA,ZINPS0,ZINPA0,KF_FS,KF_UV,KF_SCALARS,&
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
    USE LEDIR_MOD
    USE UVTVD_MOD
    USE UPDSP_MOD   ,ONLY : UPDSP
    USE UPDSPB_MOD   ,ONLY : UPDSPB
    USE MPL_MODULE      ,ONLY : MPL_BARRIER
    USE TPM_GEN         ,ONLY : LSYNC_TRANS
    USE TPM_TRANS       ,ONLY : NF_SC2, NF_SC3A, NF_SC3B
    USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX
    USE ALLOCATOR_MOD


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
    REAL(KIND=JPRBT), INTENT(IN) :: ZINPS(:), ZINPA(:)
    REAL(KIND=JPRD), INTENT(IN) :: ZINPS0(:), ZINPA0(:)

    !     LOCAL INTEGER SCALARS
    INTEGER(KIND=JPIM) :: IFC, IIFC, IDGLU, IFIRST

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    REAL(KIND=JPRB), POINTER :: POA1_L(:), POA1(:,:,:)
    REAL(KIND=JPRB), POINTER :: POA2_L(:), POA2(:,:,:)
    REAL(KIND=JPRB), POINTER :: PU(:,:,:), PV(:,:,:), PVOR(:,:,:), PDIV(:,:,:)
    REAL(KIND=JPRBT), POINTER :: ZOUT(:)
    REAL(KIND=JPRD), POINTER :: ZOUT0(:)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(LTDIR_HANDLE), INTENT(IN) :: HLTDIR
      INTEGER(KIND=C_SIZE_T) :: IALLOC_POS, IALLOC_SZ
      INTEGER(KIND=JPIM)  :: IOUT_STRIDES0, IOUT_SIZE
      INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE



    !     ------------------------------------------------------------------
    IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',0,ZHOOK_HANDLE)

    !     ------------------------------------------------------------------

    !*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
    !              --------------------------------------


    !     ------------------------------------------------------------------

    !*       2.    PREPARE WORK ARRAYS.
    !              --------------------
    CALL LEDIR_STRIDES(KF_FS,IOUT_STRIDES0=IOUT_STRIDES0,IOUT_SIZE=IOUT_SIZE,&
                       IOUT0_STRIDES0=IOUT0_STRIDES0,IOUT0_SIZE=IOUT0_SIZE)

    IALLOC_POS = 1

    IALLOC_SZ = ALIGN(2*KF_FS*(R%NTMAX+3)*D%NUMP*SIZEOF(POA1_L(1)),128)
    CALL ASSIGN_PTR(POA1_L, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ)
    CALL C_F_POINTER(C_LOC(POA1_L), POA1, (/ 2*KF_FS, R%NTMAX+3, D%NUMP /))
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    IALLOC_SZ = ALIGN(4*KF_UV*(R%NTMAX+3)*D%NUMP*SIZEOF(POA2_L(1)),128)
    CALL ASSIGN_PTR(POA2_L, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ)
    CALL C_F_POINTER(C_LOC(POA2_L), POA2, (/ 4*KF_UV, R%NTMAX+3, D%NUMP /))
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZOUT
    IALLOC_SZ = ALIGN(IOUT_SIZE*SIZEOF(ZOUT(1)),128)
    CALL ASSIGN_PTR(ZOUT, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! ZOUT0
    IALLOC_SZ = ALIGN(IOUT0_SIZE*SIZEOF(ZOUT0(1)),128)
    CALL ASSIGN_PTR(ZOUT0, GET_ALLOCATION(ALLOCATOR, HLTDIR%HOUT_AND_POA),&
        & IALLOC_POS, IALLOC_SZ)
    IALLOC_POS = IALLOC_POS + IALLOC_SZ

    ! do the legendre transform
    CALL LEDIR(ZINPS,ZINPA,ZINPS0,ZINPA0,ZOUT,ZOUT0,POA1,KF_FS)

    !$ACC DATA COPYOUT(PSPVOR,PSPDIV) IF(KF_UV > 0)
    !$ACC DATA COPYOUT(PSPSCALAR) IF(PRESENT(PSPSCALAR) .AND. KF_SCALARS > 0)
    !$ACC DATA COPYOUT(PSPSC2) IF(NF_SC2 > 0)
    !$ACC DATA COPYOUT(PSPSC3A) IF(NF_SC3A > 0)
    !$ACC DATA COPYOUT(PSPSC3B) IF(NF_SC3B > 0)

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

       ! Compute vorticity and divergence
       CALL UVTVD(KF_UV,PU,PV,PVOR,PDIV)

       ! Write back. Note, if we have UV, the contract says we *must* have VOR/DIV
       CALL UPDSPB(KF_UV,PVOR,PSPVOR,KFLDPTRUV)
       CALL UPDSPB(KF_UV,PDIV,PSPDIV,KFLDPTRUV)

    ENDIF
    !     ------------------------------------------------------------------

    !*       6.    UPDATE SPECTRAL ARRAYS.
    !              -----------------------

    ! this is on the host, so need to cp from device, Nils
    CALL UPDSP(KF_UV,KF_SCALARS,POA1,&
     & PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC)

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
