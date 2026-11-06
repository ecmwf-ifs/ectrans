#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
! (C) Copyright 1995- ECMWF.
! (C) Copyright 1995- Meteo-France.
! (C) Copyright 2022- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRLTOG_VIEW_MOD
  USE BUFFERED_ALLOCATOR_MOD, ONLY: ALLOCATION_RESERVATION_HANDLE
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: TRLTOG_VIEW, TRLTOG_HANDLE, PREPARE_TRLTOG

  TYPE TRLTOG_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HCOMBUFR_AND_COMBUFS
  END TYPE
CONTAINS
  FUNCTION PREPARE_TRLTOG(ALLOCATOR,KF_FS,KF_GP) RESULT(HTRLTOG)
    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPIB
    USE TPM_DISTR,              ONLY: D
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, RESERVE
    USE ISO_C_BINDING,          ONLY: C_SIZEOF

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP, KF_FS
    TYPE(TRLTOG_HANDLE) :: HTRLTOG

    REAL(KIND=JPRBT) :: DUMMY

    INTEGER(KIND=JPIB) :: NELEM

    NELEM = 0
    NELEM = NELEM + ALIGN(1_JPIB*KF_GP*D%NGPTOT*C_SIZEOF(DUMMY),128) ! ZCOMBUFR
    NELEM = NELEM + ALIGN(1_JPIB*KF_FS*D%NLENGTF*C_SIZEOF(DUMMY),128) !ZCOMBUFS upper obund

    HTRLTOG%HCOMBUFR_AND_COMBUFS = RESERVE(ALLOCATOR, NELEM, "HTRLTOG%HCOMBUFR_AND_COMBUFS")
  END FUNCTION PREPARE_TRLTOG

  SUBROUTINE TRLTOG_VIEW(ALLOCATOR,HTRLTOG,PREEL_REAL,KF_FS,KF_GP,KVSET, YDGP)

    !**** *trltog * - transposition of grid point data from latitudinal
    !   to column structure. This takes place between inverse
    !                 FFT and grid point calculations.
    !                 TRLTOG_VIEW is the inverse of TRGTOL

    ! Version using CUDA-aware MPI

    !     Purpose.
    !     --------


    !**   Interface.
    !     ----------
    !        *call* *trltog(...)

    !        Explicit arguments :
    !        --------------------
    !           PREEL_REAL    -  Latitudinal data ready for direct FFT (input)
    !           YDGP     -  grid points fields list    (output)
    !           KVSET    - "v-set" for each field      (input)

    !        Implicit arguments :
    !        --------------------

    !     Method.
    !     -------
    !        See documentation

    !     Externals.
    !     ----------

    !     Reference.
    !     ----------
    !        ECMWF Research Department documentation of the IFS

    !     Author.
    !     -------
    !        MPP Group *ECMWF*

    !     Modifications.
    !     --------------
    !        Original  : 95-10-01
    !        D.Dent    : 97-08-04 Reorganisation to allow NPRTRV
    !                             to differ from NPRGPEW
    !        =99-03-29= Mats Hamrud and Deborah Salmond
    !                   JUMP in FFT's changed to 1
    !                   INDEX introduced and ZCOMBUF not used for same PE
    !         01-11-23  Deborah Salmond and John Hague
    !                   LIMP_NOOLAP Option for non-overlapping message passing
    !                               and buffer packing
    !         01-12-18  Peter Towers
    !                   Improved vector performance of LTOG_PACK,LTOG_UNPACK
    !         03-0-02   G. Radnoti: Call barrier always when nproc>1
    !         08-01-01  G.Mozdzynski: cleanup
    !         09-01-02  G.Mozdzynski: use non-blocking recv and send
    !     ------------------------------------------------------------------

    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRB, JPRBT, JPIB
    USE YOMHOOK,                ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MPL_MODULE,             ONLY: MPL_WAIT, MPL_BARRIER, MPL_ABORT, MPL_RECV, MPL_SEND
    USE TPM_GEN,                ONLY: LSYNC_TRANS, NERR, LMPOFF
    USE EQ_REGIONS_MOD,         ONLY: MY_REGION_EW, MY_REGION_NS
    USE TPM_DISTR,              ONLY: D,MYSETV, MYSETW, MTAGLG,NPRCIDS,MYPROC,NPROC,NPRTRW,NPRTRV
    USE PE2SET_MOD,             ONLY: PE2SET
    USE MPL_DATA_MODULE,        ONLY: MPL_COMM_OML, JP_NON_BLOCKING_STANDARD
    USE OML_MOD,                ONLY: OML_MY_THREAD
    USE ABORT_TRANS_MOD,        ONLY: ABORT_TRANS
#ifdef USE_RAW_MPI
    USE MPI_F08,                ONLY: MPI_COMM, MPI_REQUEST, MPI_REAL4, MPI_REAL8
    ! Missing: MPI_ISEND, MPI_IRECV on purpose due to cray-mpi bug (see https://github.com/ecmwf-ifs/ectrans/pull/157)
#endif
    USE TPM_STATS,              ONLY: GSTATS => GSTATS_NVTX
    USE TPM_TRANS,              ONLY: LDIVGP, LSCDERS, LUVDER, LVORGP, NPROMA
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, ASSIGN_PTR, GET_ALLOCATION
    USE ISO_C_BINDING,          ONLY: C_SIZEOF
    USE OPENACC_EXT,            ONLY: EXT_ACC_ARR_DESC, EXT_ACC_PASS, EXT_ACC_CREATE, &
      &                               EXT_ACC_DELETE
#ifdef ACCGPU
    USE OPENACC,                ONLY: ACC_HANDLE_KIND
#endif
    USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : GRID_VIEW

    IMPLICIT NONE

    REAL(KIND=JPRBT),  INTENT(INOUT), POINTER  :: PREEL_REAL(:)
    INTEGER(KIND=JPIM),INTENT(IN)  :: KF_FS,KF_GP
    INTEGER(KIND=JPIM), INTENT(IN) :: KVSET(KF_GP)
    TYPE(GRID_VIEW):: YDGP(:)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(TRLTOG_HANDLE) :: HTRLTOG

    ! LOCAL VARIABLES

    REAL(KIND=JPRBT), POINTER :: ZCOMBUFS(:),ZCOMBUFR(:)

    LOGICAL :: LLOCAL_CONTRIBUTION
    INTEGER(KIND=JPIB) :: ISENDTOT (NPROC)
    INTEGER(KIND=JPIB) :: IRECVTOT (NPROC)
    INTEGER(KIND=JPIM) :: ISENDTOT_MPI(NPROC)
    INTEGER(KIND=JPIM) :: IRECVTOT_MPI(NPROC)
    INTEGER(KIND=JPIM) :: IREQ     (NPROC*2)
    INTEGER(KIND=JPIM) :: IRECV_TO_PROC(NPROC)
    INTEGER(KIND=JPIM) :: ISEND_TO_PROC(NPROC)

    INTEGER(KIND=JPIM) :: JFLD, J, JI, JGL, JK, JL, IFLDS, JROC, INR, INS
    INTEGER(KIND=JPIM) :: IFIRSTLAT, ILASTLAT, IFLD, IGL, IGLL,&
                 &ISETA, ISETB, ISETV, ISEND, IRECV, ISETW, IPROC, &
                 &IR, ILOCAL_LAT, ISEND_COUNTS, IRECV_COUNTS, IERROR, II, ILEN, &
                 &JBLK, ILAT_STRIP
    INTEGER(KIND=JPIB) :: IPOS

    ! Contains FIELD, PARS, LEVS
    INTEGER(KIND=JPIM) :: IUVPAR,IGP2PAR,IGP3ALEV,IGP3APAR,IGP3BLEV,IGP3BPAR

    INTEGER(KIND=JPIB) :: IIN_TO_SEND_BUFR(D%NLENGTF,2)
    INTEGER(KIND=JPIM) :: IIN_TO_SEND_BUFR_OFFSET(NPROC), IIN_TO_SEND_BUFR_V
    INTEGER(KIND=JPIM) :: IRECV_FIELD_COUNT(NPRTRV),IRECV_FIELD_COUNT_V
    INTEGER(KIND=JPIM) :: IRECV_WSET_SIZE(NPRTRW),IRECV_WSET_SIZE_V
    INTEGER(KIND=JPIM) :: IRECV_WSET_OFFSET(NPRTRW+1), IRECV_WSET_OFFSET_V
    INTEGER(KIND=JPIB), ALLOCATABLE :: ICOMBUFS_OFFSET(:),ICOMBUFR_OFFSET(:), IFLDA(:,:)
    INTEGER(KIND=JPIB) :: ICOMBUFS_OFFSET_V, ICOMBUFR_OFFSET_V
    
    INTEGER(KIND=JPIM) :: J3

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
        
#ifdef USE_RAW_MPI
    TYPE(MPI_COMM) :: LOCAL_COMM
    TYPE(MPI_REQUEST) :: IREQUEST(NPROC*2)
#else
    INTEGER(KIND=JPIM) :: IREQUEST(NPROC*2)
#endif

#ifdef PARKINDTRANS_SINGLE
#define TRLTOG_DTYPE MPI_REAL4
#else
#define TRLTOG_DTYPE MPI_REAL8
#endif
#ifdef USE_RAW_MPI
    IF(.NOT. LMPOFF) THEN
      LOCAL_COMM%MPI_VAL = MPL_COMM_OML( OML_MY_THREAD() )
    ENDIF
#endif
    !     ------------------------------------------------------------------

    !*       0.    Some initializations
    !              --------------------
    IF (LHOOK) CALL DR_HOOK('TRLTOG_VIEW',0,ZHOOK_HANDLE)

    ! We first get the decomposition individually
        
    CALL GSTATS(1806,0)

    ! Prepare receiver arrays
    ! find number of fields on a certain V-set
    IF(NPRTRV == 1) THEN
      ! This is needed because KVSET(JFLD) == -1 if there is only one V-set
      IRECV_FIELD_COUNT(1) = KF_GP
    ELSE
      IRECV_FIELD_COUNT(:) = 0
      DO JFLD=1,KF_GP
        IRECV_FIELD_COUNT(KVSET(JFLD)) = IRECV_FIELD_COUNT(KVSET(JFLD)) + 1
      ENDDO
    ENDIF
    ! find number of grid-points on a certain W-set that overlap with myself
    IRECV_WSET_SIZE(:) = 0
    DO ILOCAL_LAT=D%NFRSTLAT(MY_REGION_NS),D%NLSTLAT(MY_REGION_NS)
      ILAT_STRIP = ILOCAL_LAT-D%NFRSTLAT(MY_REGION_NS)+D%NPTRFLOFF+1
      IRECV_WSET_SIZE(D%NPROCL(ILOCAL_LAT)) = &
          & IRECV_WSET_SIZE(D%NPROCL(ILOCAL_LAT))+D%NONL(ILAT_STRIP,MY_REGION_EW)
    ENDDO
    ! sum up offsets
    IRECV_WSET_OFFSET(1) = 0
    DO JROC=1,NPRTRW
      IRECV_WSET_OFFSET(JROC+1)=IRECV_WSET_OFFSET(JROC)+IRECV_WSET_SIZE(JROC)
    ENDDO
    DO JROC=1,NPROC
      CALL PE2SET(JROC,ISETA,ISETB,ISETW,ISETV)
      ! total recv size is # points per field * # fields
      IRECVTOT(JROC) = 1_JPIB*IRECV_WSET_SIZE(ISETW)*IRECV_FIELD_COUNT(ISETV)
    ENDDO

    ! Prepare sender arrays
    IIN_TO_SEND_BUFR_OFFSET(1) = 0
    DO JROC=1,NPROC
      ! Get new offset to my current KINDEX entry
      IF (JROC > 1 .AND. KF_FS > 0) THEN
        IIN_TO_SEND_BUFR_OFFSET(JROC) = IIN_TO_SEND_BUFR_OFFSET(JROC-1)+ISENDTOT(JROC-1)/KF_FS
      ELSEIF (JROC > 1) THEN
        IIN_TO_SEND_BUFR_OFFSET(JROC) = IIN_TO_SEND_BUFR_OFFSET(JROC-1)
      ENDIF

      CALL PE2SET(JROC,ISETA,ISETB,ISETW,ISETV)

      ! MAX(Index of first fourier latitude for this W set, first latitude of a senders A set)
      ! i.e. we find the overlap between what we have on sender side (others A set) and the receiver
      ! (me, the W-set). Ideally those conincide, at least mostly.
      IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(ISETA))
      ! MIN(Index of last fourier latitude for this W set, last latitude of a senders A set)
      ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(ISETA))

      IPOS = 0
      DO JGL=IFIRSTLAT,ILASTLAT
        ! get from "actual" latitude to the latitude strip offset
        IGL  = JGL-D%NFRSTLAT(ISETA)+D%NPTRFRSTLAT(ISETA)
        ! get from "actual" latitude to the latitude offset
        IGLL = JGL-D%NPTRLS(MYSETW)+1
        DO JL=1,D%NONL(IGL,ISETB)
          IPOS = IPOS+1
          ! offset to first layer of this gridpoint
          IIN_TO_SEND_BUFR(IIN_TO_SEND_BUFR_OFFSET(JROC)+IPOS,1) = &
              & 1_JPIB*KF_FS*D%NSTAGTF(IGLL)+(D%NSTA(IGL,ISETB)-1)+(JL-1)
          ! distance between two layers of this gridpoint
          IIN_TO_SEND_BUFR(IIN_TO_SEND_BUFR_OFFSET(JROC)+IPOS,2) = &
              & D%NSTAGTF(IGLL+1)-D%NSTAGTF(IGLL)
        ENDDO
      ENDDO
      !we always receive the full fourier space
      ISENDTOT(JROC) = IPOS*KF_FS
    ENDDO
    LLOCAL_CONTRIBUTION = ISENDTOT(MYPROC) > 0
    
#ifdef OMPGPU    
    !$OMP TARGET DATA MAP(TO:IIN_TO_SEND_BUFR) MAP(PRESENT,ALLOC:PREEL_REAL) IF(KF_FS > 0)
#endif
#ifdef ACCGPU    
    ! Present until self contribution and packing are done
    !$ACC DATA COPYIN(IIN_TO_SEND_BUFR) PRESENT(PREEL_REAL) IF(KF_FS > 0) ASYNC(1)
#endif

    CALL GSTATS(1806,1)
    
    ! Figure out processes that send or recv something
    ISEND_COUNTS   = 0
    IRECV_COUNTS   = 0
    DO JROC=1,NPROC
      IF( JROC /= MYPROC) THEN
        IF(IRECVTOT(JROC) > 0) THEN
          ! I have to recv something, so let me store that
          IRECV_COUNTS = IRECV_COUNTS + 1
          IRECV_TO_PROC(IRECV_COUNTS)=JROC
        ENDIF
        IF(ISENDTOT(JROC) > 0) THEN
          ! I have to send something, so let me store that
          ISEND_COUNTS = ISEND_COUNTS+1
          ISEND_TO_PROC(ISEND_COUNTS)=JROC
        ENDIF
      ENDIF
    ENDDO

    ! ... build this data structure now during the MPI communication
    ! Allocate this buffer now. Add 1 for self contribution
    ALLOCATE(IFLDA(KF_GP,1+IRECV_COUNTS))

    ! Copy local contribution
    IF(LLOCAL_CONTRIBUTION) THEN
      ! I have to send something to myself...

      ! Input is KF_GP fields. We find the resulting KF_FS fields.
      IFLDS = 0
      DO JFLD=1,KF_GP
        IF(KVSET(JFLD) == MYSETV .OR. KVSET(JFLD) == -1) THEN
          IFLDS = IFLDS+1          
          IFLDA(IFLDS,1) = JFLD          
        ENDIF
      ENDDO
    ENDIF

    DO INR=1,IRECV_COUNTS
      IRECV=IRECV_TO_PROC(INR)
      CALL PE2SET(IRECV,ISETA,ISETB,ISETW,ISETV)
      IFLDS = 0
      DO JFLD=1,KF_GP
        IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1 ) THEN
          IFLDS = IFLDS+1          
          IFLDA(IFLDS,1+INR)=JFLD          
        ENDIF
      ENDDO
    ENDDO
   
#ifdef OMPGPU
    !$OMP TARGET DATA MAP(TO:IFLDA)
#endif
#ifdef ACCGPU
    !$ACC DATA COPYIN(IFLDA,YDGP) ASYNC(1)
#endif

    ! Copy local contribution
    IF(LLOCAL_CONTRIBUTION) THEN

      CALL GSTATS(1604,0)

      IRECV_WSET_OFFSET_V = IRECV_WSET_OFFSET(MYSETW)
      IRECV_WSET_SIZE_V = IRECV_WSET_SIZE(MYSETW)
      IIN_TO_SEND_BUFR_V = IIN_TO_SEND_BUFR_OFFSET(MYPROC)
      
#ifdef OMPGPU
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
      !$OMP& PRIVATE(JK,JBLK,IFLD,IPOS) &
      !$OMP& SHARED(KF_FS,IRECV_WSET_SIZE_V,NPROMA,IRECV_WSET_OFFSET_V,IFLDA,IIN_TO_SEND_BUFR_V, &
      !$OMP&        IIN_TO_SEND_BUFR,PREEL_REAL,YDGP) &
      !$OMP& MAP(TO:KF_FS,IRECV_WSET_SIZE_V,NPROMA,IRECV_WSET_OFFSET_V,IIN_TO_SEND_BUFR_V)
#endif
#ifdef ACCGPU
      !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(JK,JBLK,IFLD,IPOS) &
      !$ACC&         FIRSTPRIVATE(KF_FS,IRECV_WSET_SIZE_V,IRECV_WSET_OFFSET_V, &
      !$ACC&         IIN_TO_SEND_BUFR_V,NPROMA) ASYNC(1)
#endif
      DO JFLD=1,KF_FS
        DO JL=1,IRECV_WSET_SIZE_V
          JK = MOD(IRECV_WSET_OFFSET_V+JL-1,NPROMA)+1
          JBLK = (IRECV_WSET_OFFSET_V+JL-1)/NPROMA+1
          IFLD = IFLDA(JFLD,1)
          IPOS = IIN_TO_SEND_BUFR(IIN_TO_SEND_BUFR_V+JL,1)+ &
              & (JFLD-1)*IIN_TO_SEND_BUFR(IIN_TO_SEND_BUFR_V+JL,2)+1
          YDGP(IFLD)%P(JK,JBLK) = PREEL_REAL(IPOS)
        ENDDO
      ENDDO
    
      CALL GSTATS(1604,1)
    ENDIF

    ALLOCATE(ICOMBUFS_OFFSET(ISEND_COUNTS+1))
    ICOMBUFS_OFFSET(1) = 0
    DO JROC=1,ISEND_COUNTS
      ICOMBUFS_OFFSET(JROC+1) = ICOMBUFS_OFFSET(JROC) + ISENDTOT(ISEND_TO_PROC(JROC))
    ENDDO
    ALLOCATE(ICOMBUFR_OFFSET(IRECV_COUNTS+1))
    ICOMBUFR_OFFSET(1) = 0
    DO JROC=1,IRECV_COUNTS
      ICOMBUFR_OFFSET(JROC+1) = ICOMBUFR_OFFSET(JROC) + IRECVTOT(IRECV_TO_PROC(JROC))
    ENDDO

    IF (IRECV_COUNTS > 0) THEN
      CALL ASSIGN_PTR(ZCOMBUFR, GET_ALLOCATION(ALLOCATOR, HTRLTOG%HCOMBUFR_AND_COMBUFS),&
          & 1_JPIB, ICOMBUFR_OFFSET(IRECV_COUNTS+1)*C_SIZEOF(ZCOMBUFR(1)))
    ENDIF
    IF (ISEND_COUNTS > 0) THEN
      CALL ASSIGN_PTR(ZCOMBUFS, GET_ALLOCATION(ALLOCATOR, HTRLTOG%HCOMBUFR_AND_COMBUFS),&
          & ALIGN(1_JPIB*KF_GP*D%NGPTOT*C_SIZEOF(ZCOMBUFR(1)),128)+1, &
          & ICOMBUFS_OFFSET(ISEND_COUNTS+1)*C_SIZEOF(ZCOMBUFS(1)))
    ENDIF

#ifdef OMPGPU
    !$OMP TARGET DATA MAP(PRESENT,ALLOC:ZCOMBUFS) IF(ISEND_COUNTS > 0)
#endif
#ifdef ACCGPU
    !$ACC DATA PRESENT(ZCOMBUFS) IF(ISEND_COUNTS > 0) ASYNC(1)
#endif
    CALL GSTATS(1605,0)
    DO INS=1,ISEND_COUNTS
      IPROC = ISEND_TO_PROC(INS)
      ILEN = ISENDTOT(IPROC)/KF_FS
      IIN_TO_SEND_BUFR_V = IIN_TO_SEND_BUFR_OFFSET(IPROC)
      ICOMBUFS_OFFSET_V = ICOMBUFS_OFFSET(INS)
#ifdef OMPGPU
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
      !$OMP& PRIVATE(IPOS) &
      !$OMP& SHARED(KF_FS,ILEN,IIN_TO_SEND_BUFR_V,IIN_TO_SEND_BUFR,PREEL_REAL,ICOMBUFS_OFFSET_V, &
      !$OMP&        ZCOMBUFS) &
      !$OMP& MAP(TO:KF_FS,ILEN,IIN_TO_SEND_BUFR_V,ICOMBUFS_OFFSET_V)
#endif
#ifdef ACCGPU
      !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IPOS) FIRSTPRIVATE(KF_FS,ILEN,IIN_TO_SEND_BUFR_V, &
      !$ACC&              ICOMBUFS_OFFSET_V) COLLAPSE(2) ASYNC(1)
#endif
      DO JFLD=1,KF_FS
        DO JL=1,ILEN
          IPOS = IIN_TO_SEND_BUFR(IIN_TO_SEND_BUFR_V+JL,1)+ &
              & (JFLD-1)*IIN_TO_SEND_BUFR(IIN_TO_SEND_BUFR_V+JL,2)+1
          ZCOMBUFS(ICOMBUFS_OFFSET_V+(JFLD-1)*ILEN+JL) = PREEL_REAL(IPOS)
        ENDDO
      ENDDO
    ENDDO
    CALL GSTATS(1605,1)
#ifdef OMPGPU
    !$OMP END TARGET DATA ! ZCOMBUFS
#endif
#ifdef ACCGPU
    !$ACC END DATA ! ZCOMBUFS

    !$ACC WAIT(1)
#endif

    CALL GSTATS(805,0)

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(421,0)

    IR=0
    !...Receive loop.........................................................
#ifdef USE_GPU_AWARE_MPI
#ifdef OMPGPU
    !$OMP TARGET DATA USE_DEVICE_PTR(ZCOMBUFS,ZCOMBUFR)
#endif
#ifdef ACCGPU
    !$ACC HOST_DATA USE_DEVICE(ZCOMBUFS,ZCOMBUFR)
#endif
#else
#ifdef OMPGPU
    !$OMP TARGET UPDATE FROM(ZCOMBUFS) IF(ISEND_COUNTS > 0)
#endif
#ifdef ACCGPU
    !! this is safe-but-slow fallback for running without GPU-aware MPI
    !$ACC UPDATE HOST(ZCOMBUFS) IF(ISEND_COUNTS > 0)
#endif
#endif

    ! Skip the own contribution because this is ok to overflow
    ISENDTOT(MYPROC) = 0
    IRECVTOT(MYPROC) = 0

    ISENDTOT_MPI = ISENDTOT
    IRECVTOT_MPI = IRECVTOT
    IF (ANY(ISENDTOT_MPI /= ISENDTOT)) &
      & CALL MPL_ABORT("Overflow in trltog")
    IF (ANY(IRECVTOT_MPI /= IRECVTOT)) &
      & CALL MPL_ABORT("Overflow in trltog")

    DO INR=1,IRECV_COUNTS
      IR=IR+1
      IRECV=IRECV_TO_PROC(INR)
#ifdef USE_RAW_MPI
      CALL MPI_IRECV(ZCOMBUFR(ICOMBUFR_OFFSET(INR)+1:ICOMBUFR_OFFSET(INR+1)), &
        & IRECVTOT_MPI(IRECV), &
        & TRLTOG_DTYPE,NPRCIDS(IRECV)-1, &
        & MTAGLG, LOCAL_COMM, IREQUEST(IR), &
        & IERROR )
      IREQ(IR) = IREQUEST(IR)%MPI_VAL
#else
      CALL MPL_RECV(ZCOMBUFR(ICOMBUFR_OFFSET(INR)+1:ICOMBUFR_OFFSET(INR+1)), &
        &           KSOURCE=NPRCIDS(IRECV), KTAG=MTAGLG, KMP_TYPE=JP_NON_BLOCKING_STANDARD, &
        &           KREQUEST=IREQUEST(IR))
      IREQ(IR) = IREQUEST(IR)
#endif
    ENDDO

    !...Send loop.........................................................
    DO INS=1,ISEND_COUNTS
      IR=IR+1
      ISEND=ISEND_TO_PROC(INS)
#ifdef USE_RAW_MPI
      CALL MPI_ISEND(ZCOMBUFS(ICOMBUFS_OFFSET(INS)+1:ICOMBUFS_OFFSET(INS+1)),ISENDTOT_MPI(ISEND), &
        & TRLTOG_DTYPE, NPRCIDS(ISEND)-1,MTAGLG,LOCAL_COMM,IREQUEST(IR),IERROR)
      IREQ(IR) = IREQUEST(IR)%MPI_VAL
#else
      CALL MPL_SEND(ZCOMBUFS(ICOMBUFS_OFFSET(INS)+1:ICOMBUFS_OFFSET(INS+1)), &
        &           KDEST=NPRCIDS(ISEND), KTAG=MTAGLG, KMP_TYPE=JP_NON_BLOCKING_STANDARD, &
        &           KREQUEST=IREQUEST(IR))
      IREQ(IR) = IREQUEST(IR)
#endif
    ENDDO

    IF(IR > 0) THEN
      CALL MPL_WAIT(KREQUEST=IREQ(1:IR), &
      & CDSTRING='TRLTOG_VIEW: WAIT FOR SENDS AND RECEIVES')
    ENDIF

#ifdef USE_GPU_AWARE_MPI
#ifdef ACCGPU
    !$ACC END HOST_DATA ! ZCOMBUFS, ZCOMBUFR
#endif
#ifdef OMPGPU
    !$OMP END TARGET DATA ! ZCOMBUFS, ZCOMBUFR
#endif
#else
#ifdef OMPGPU
#endif
    !! this is safe-but-slow fallback for running without GPU-aware MPI
#ifdef OMPGPU
    !$OMP TARGET UPDATE TO(ZCOMBUFR) IF(IRECV_COUNTS > 0)
#endif
#ifdef ACCGPU
    !$ACC UPDATE DEVICE(ZCOMBUFR) IF(IRECV_COUNTS > 0)
#endif
#endif

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(441,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(441,1)
    ENDIF
    CALL GSTATS(421,1)

#ifdef OMPGPU
    !$OMP TARGET DATA MAP(PRESENT,ALLOC:ZCOMBUFR) IF(IRECV_COUNTS > 0)
#endif
#ifdef ACCGPU
    !$ACC DATA PRESENT(ZCOMBUFR) IF(IRECV_COUNTS > 0) ASYNC(1)
#endif
    CALL GSTATS(805,1)

    !  Unpack loop.........................................................

    CALL GSTATS(1606,0)
    DO INR=1,IRECV_COUNTS
      IRECV=IRECV_TO_PROC(INR)
      CALL PE2SET(IRECV,ISETA,ISETB,ISETW,ISETV)

      IRECV_FIELD_COUNT_V = IRECV_FIELD_COUNT(ISETV)
      ICOMBUFR_OFFSET_V = ICOMBUFR_OFFSET(INR)

      IRECV_WSET_OFFSET_V = IRECV_WSET_OFFSET(ISETW)
      IRECV_WSET_SIZE_V = IRECV_WSET_SIZE(ISETW)
      
#ifdef OMPGPU
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
      !$OMP& PRIVATE(JK,JBLK,IFLD,JI) &
      !$OMP& SHARED(IRECV_FIELD_COUNT_V,IRECV_WSET_SIZE_V,NPROMA,IRECV_WSET_OFFSET_V,IFLDA, &
      !$OMP&        ICOMBUFR_OFFSET_V,ZCOMBUFR,INR,YDGP) &
      !$OMP& MAP(TO:IRECV_FIELD_COUNT_V,IRECV_WSET_SIZE_V,NPROMA,IRECV_WSET_OFFSET_V, &
      !$OMP&     ICOMBUFR_OFFSET_V)
#endif
#ifdef ACCGPU
      !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(JK,JBLK,IFLD,JI) &
      !$ACC&              FIRSTPRIVATE(INR,IRECV_FIELD_COUNT_V,IRECV_WSET_SIZE_V,&
      !$ACC&              IRECV_WSET_OFFSET_V,NPROMA,ICOMBUFR_OFFSET_V) ASYNC(1)
#endif
      DO JFLD=1,IRECV_FIELD_COUNT_V
        DO JL=1,IRECV_WSET_SIZE_V
          JK = MOD(IRECV_WSET_OFFSET_V+JL-1,NPROMA)+1
          JBLK = (IRECV_WSET_OFFSET_V+JL-1)/NPROMA+1
          IFLD=IFLDA(JFLD,1+INR)
          JI = ICOMBUFR_OFFSET_V+(JFLD-1)*IRECV_WSET_SIZE_V+JL
          YDGP(IFLD)%P(JK,JBLK) = ZCOMBUFR(JI)
        ENDDO
      ENDDO      
    
    ENDDO
#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC WAIT(1)
#endif

#ifdef OMPGPU
    !$OMP END TARGET DATA ! ZCOMBUFR
#endif
#ifdef ACCGPU
    !$ACC END DATA ! ZCOMBUFR
#endif
    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(422,0)
#ifdef OMPGPU
    !$OMP END TARGET DATA ! IFLDA
    !$OMP END TARGET DATA ! PREEL_REAL    
#endif
#ifdef ACCGPU
    !$ACC END DATA ! IFLDA
    !$ACC END DATA ! PREEL_REAL    
#endif
    
    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(442,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(442,1)
    ENDIF
    CALL GSTATS(422,1)

    CALL GSTATS(1606,1)

    ! Free this now
    DEALLOCATE(IFLDA)

    IF (LHOOK) CALL DR_HOOK('TRLTOG_VIEW',1,ZHOOK_HANDLE)
  END SUBROUTINE TRLTOG_VIEW
END MODULE TRLTOG_VIEW_MOD

