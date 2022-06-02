! (C) Copyright 1995- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRGTOL_MOD
  CONTAINS
  SUBROUTINE TRGTOL_CUDAAWARE(PREEL_REAL,KF_FS,KF_GP,KVSET,KPTRGP,&
   &PGP,PGPUV,PGP3A,PGP3B,PGP2)

  !**** *TRGTOL * - transposition of grid point data from column
  !                 structure to latitudinal. Reorganize data between
  !                 grid point calculations and direct Fourier Transform

  ! Version using CUDA-aware MPI

  !     Purpose.
  !     --------


  !**   Interface.
  !     ----------
  !        *call* *trgtol(...)

  !        Explicit arguments :
  !        --------------------
  !           PREEL_REAL    -  Latitudinal data ready for direct FFT (output)
  !           PGP    -  Blocked grid point data    (input)

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
  !        Original: 95-10-01
  !        D.Dent  : 97-08-04   Reorganisation to allow
  !                             NPRTRV to differ from NPRGPEW
  !                : 98-06-17   add mailbox control logic (from TRLTOM)
  !        =99-03-29= Mats Hamrud and Deborah Salmond
  !                   JUMP in FFT's changed to 1
  !                   KINDEX introduced and ZCOMBUF not used for same PE
  !         01-11-23  Deborah Salmond and John Hague
  !                    LIMP_NOOLAP Option for non-overlapping message passing
  !                    and buffer packing
  !         01-12-18  Peter Towers
  !                   Improved vector performance of GTOL_PACK,GTOL_UNPACK
  !         03-04-02  G. Radnoti: call barrier always when nproc>1
  !         08-01-01  G.Mozdzynski: cleanup
  !         09-01-02  G.Mozdzynski: use non-blocking recv and send
  !     ------------------------------------------------------------------



  USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRB ,  JPRBT
  USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK,  JPHOOK
  USE MPL_MODULE      ,ONLY : MPL_WAIT, MPL_BARRIER
  USE TPM_GEN         ,ONLY : LSYNC_TRANS
  USE EQ_REGIONS_MOD  ,ONLY : MY_REGION_EW, MY_REGION_NS
  USE TPM_DISTR       ,ONLY : D,MYSETV, MYSETW, MTAGLG,NPRCIDS,MYPROC,NPROC,NPRTRW,NPRTRV
  USE PE2SET_MOD      ,ONLY : PE2SET
  USE MPL_DATA_MODULE ,ONLY : MPL_COMM_OML
  USE OML_MOD         ,ONLY : OML_MY_THREAD
  USE MPI
  USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

  USE TPM_TRANS       ,ONLY : NPROMA

  IMPLICIT NONE

  REAL(KIND=JPRBT),INTENT(OUT)   :: PREEL_REAL(:)
  INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(:)
  INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
  INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
  REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP(:,:,:)
  REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGPUV(:,:,:,:)
  REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3A(:,:,:,:)
  REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3B(:,:,:,:)
  REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP2(:,:,:)

  ! LOCAL VARIABLES

  !     LOCAL INTEGER SCALARS
  REAL(KIND=JPRBT),ALLOCATABLE :: ZCOMBUFS(:),ZCOMBUFR(:)

  INTEGER(KIND=JPIM) :: ISENDTOT (NPROC)
  INTEGER(KIND=JPIM) :: IRECVTOT (NPROC)
  INTEGER(KIND=JPIM) :: IREQ     (NPROC*2)
  INTEGER(KIND=JPIM) :: IRECV_TO_PROC(NPROC)
  INTEGER(KIND=JPIM) :: ISEND_TO_PROC(NPROC)

  INTEGER(KIND=JPIM) :: IFIRSTLAT, IGL, IGLL, ILAST,&
               &ILASTLAT, ILEN, JROC, IPOS, ISETA, &
               &ISETB, IRECV, &
               &ISETV, ISEND, JBLK, JFLD, &
               &JGL, JI, JK, JL, ISETW,  IFLD, &
               &II,IBUFLENR,IRECV_COUNTS, IPROC,IFLDS, &
               &ISEND_COUNTS,INS,INR,IR, JKL, PBOUND, IERROR, ILOCAL_LAT
  INTEGER(KIND=JPIM) :: KF, KGL, KI

  INTEGER(KIND=JPIM) :: IOFF, ILAT_STRIP
  INTEGER(KIND=JPIM) :: IRECV_BUFR_TO_OUT(D%NLENGTF,2),IRECV_BUFR_TO_OUT_OFFSET(NPROC), IRECV_BUFR_TO_OUT_V
  INTEGER(KIND=JPIM) :: ISEND_FIELD_COUNT(NPRTRV),ISEND_FIELD_COUNT_V
  INTEGER(KIND=JPIM) :: ISEND_WSET_SIZE(NPRTRW),ISEND_WSET_SIZE_V
  INTEGER(KIND=JPIM) :: ISEND_WSET_OFFSET(NPRTRW+1), ISEND_WSET_OFFSET_V
  INTEGER(KIND=JPIM), ALLOCATABLE :: ICOMBUFS_OFFSET(:),ICOMBUFR_OFFSET(:)
  INTEGER(KIND=JPIM) :: ICOMBUFS_OFFSET_V, ICOMBUFR_OFFSET_V
  INTEGER(KIND=JPIM) :: IFLDA(KF_GP)

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR

  INTEGER(JPIM), PARAMETER :: PGP_INDICES_UV = 1
  INTEGER(JPIM), PARAMETER :: PGP_INDICES_GP2 = 2
  INTEGER(JPIM), PARAMETER :: PGP_INDICES_GP3A = 3
  INTEGER(JPIM), PARAMETER :: PGP_INDICES_GP3B = 4
  INTEGER(JPIM), PARAMETER :: PGP_INDICES_END = 5
  INTEGER(JPIM) :: PGP_INDICES(PGP_INDICES_END)

  #ifdef PARKINDTRANS_SINGLE
  #define TRGTOL_DTYPE MPI_REAL
  #else
  #define TRGTOL_DTYPE MPI_DOUBLE_PRECISION
  #endif

  !     ------------------------------------------------------------------

  !*       0.    Some initializations
  !              --------------------

  IF (LHOOK) CALL DR_HOOK('TRGTOL_CUDAAWARE',0,ZHOOK_HANDLE)

  CALL GSTATS(1805,0)
  IOFF=1
  PGP_INDICES(PGP_INDICES_UV) = IOFF
  IF (PRESENT(PGPUV)) IOFF=IOFF+UBOUND(PGPUV,2)*2
  PGP_INDICES(PGP_INDICES_GP2) = IOFF
  IF (PRESENT(PGP2)) IOFF=IOFF+UBOUND(PGP2,2)
  PGP_INDICES(PGP_INDICES_GP3A) = IOFF
  IF (PRESENT(PGP3A)) IOFF=IOFF+UBOUND(PGP3A,2)*UBOUND(PGP3A,3)
  PGP_INDICES(PGP_INDICES_GP3B) = IOFF
  IF (PRESENT(PGP3B)) IOFF=IOFF+UBOUND(PGP3B,2)*UBOUND(PGP3B,3)
  PGP_INDICES(PGP_INDICES_END) = IOFF

  ! Prepare sender arrays
  ! find number of fields on a certain V-set
  IF(NPRTRV == 1) THEN
    ! This is needed because KVSET(JFLD) == -1 if there is only one V-set
    ISEND_FIELD_COUNT(1) = KF_GP
  ELSE
    ISEND_FIELD_COUNT(:) = 0
    DO JFLD=1,KF_GP
      ISEND_FIELD_COUNT(KVSET(JFLD)) = ISEND_FIELD_COUNT(KVSET(JFLD)) + 1
    ENDDO
  ENDIF
  ! find number of grid-points on a certain W-set that overlap with myself
  ISEND_WSET_SIZE(:) = 0
  DO ILOCAL_LAT=D%NFRSTLAT(MY_REGION_NS),D%NLSTLAT(MY_REGION_NS)
    ILAT_STRIP = ILOCAL_LAT-D%NFRSTLAT(MY_REGION_NS)+D%NPTRFLOFF+1
    ISEND_WSET_SIZE(D%NPROCL(ILOCAL_LAT)) = &
        & ISEND_WSET_SIZE(D%NPROCL(ILOCAL_LAT))+D%NONL(ILAT_STRIP,MY_REGION_EW)
  ENDDO
  ! sum up offsets
  ISEND_WSET_OFFSET(1) = 0
  DO JROC=1,NPRTRW
    ISEND_WSET_OFFSET(JROC+1)=ISEND_WSET_OFFSET(JROC)+ISEND_WSET_SIZE(JROC)
  ENDDO
  DO JROC=1,NPROC
    CALL PE2SET(JROC,ISETA,ISETB,ISETW,ISETV)
    ! total send size is # points per field * # fields
    ISENDTOT(JROC) = ISEND_WSET_SIZE(ISETW)*ISEND_FIELD_COUNT(ISETV)
  ENDDO

  ! Prepare receiver arrays
  IRECV_BUFR_TO_OUT_OFFSET(:) = 0
  DO JROC=1,NPROC
    ! Get new offset to my current KINDEX entry
    IF (JROC > 1 .AND. KF_FS > 0) THEN
      IRECV_BUFR_TO_OUT_OFFSET(JROC) = IRECV_BUFR_TO_OUT_OFFSET(JROC-1)+IRECVTOT(JROC-1)/KF_FS
    ELSEIF (JROC > 1) THEN
      IRECV_BUFR_TO_OUT_OFFSET(JROC) = IRECV_BUFR_TO_OUT_OFFSET(JROC-1)
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
        IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_OFFSET(JROC)+IPOS,1) = &
            & KF_FS*D%NSTAGTF(IGLL)+(D%NSTA(IGL,ISETB)-1)+(JL-1)
        ! distance between two layers of this gridpoint
        IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_OFFSET(JROC)+IPOS,2) = &
            & D%NSTAGTF(IGLL+1)-D%NSTAGTF(IGLL)
      ENDDO
    ENDDO
    !we always receive the full fourier space
    IRECVTOT(JROC) = IPOS*KF_FS
  ENDDO

  !$ACC DATA COPYIN(IRECV_BUFR_TO_OUT,PGP_INDICES) PRESENT(PREEL_REAL) ASYNC(1)

  CALL GSTATS(1805,1)

  ! Put data on device for copyin
  IF (LSYNC_TRANS) THEN
    !$ACC WAIT(1)
    CALL GSTATS(430,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(430,1)
  ENDIF
  CALL GSTATS(412,0)
  !$ACC DATA IF(PRESENT(PGP))   COPYIN(PGP) ASYNC(1)
  !$ACC DATA IF(PRESENT(PGPUV)) COPYIN(PGPUV) ASYNC(1)
  !$ACC DATA IF(PRESENT(PGP2))  COPYIN(PGP2) ASYNC(1)
  !$ACC DATA IF(PRESENT(PGP3A)) COPYIN(PGP3A) ASYNC(1)
  !$ACC DATA IF(PRESENT(PGP3B)) COPYIN(PGP3B) ASYNC(1)
  IF (LSYNC_TRANS) THEN
    !$ACC WAIT(1)
    CALL GSTATS(432,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(432,1)
  ENDIF
  CALL GSTATS(412,1)

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

  ! Do this with "enter data" syntax because we are in the PGP data clause
  IF (ISEND_COUNTS > 0) ALLOCATE(ZCOMBUFS(ICOMBUFS_OFFSET(ISEND_COUNTS+1)))
  !$ACC ENTER DATA IF(ISEND_COUNTS > 0) CREATE(ZCOMBUFS) ASYNC(1)

  !....Pack loop.........................................................
  !$ACC DATA IF(ISEND_COUNTS > 0) PRESENT(ZCOMBUFS) ASYNC(1)

  CALL GSTATS(1602,0)
  DO INS=1,ISEND_COUNTS
    ISEND=ISEND_TO_PROC(INS)
    CALL PE2SET(ISEND,ISETA,ISETB,ISETW,ISETV)

    ISEND_FIELD_COUNT_V = ISEND_FIELD_COUNT(ISETV)
    ICOMBUFS_OFFSET_V = ICOMBUFS_OFFSET(INS)

    IFLDS = 0
    DO JFLD=1,KF_GP
      IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1 ) THEN
        IFLDS = IFLDS+1
        IF(PRESENT(KPTRGP)) THEN
          IFLDA(IFLDS)=KPTRGP(JFLD)
        ELSE
          IFLDA(IFLDS)=JFLD
        ENDIF
      ENDIF
    ENDDO

    !$ACC DATA COPYIN(IFLDA(1:ISEND_FIELD_COUNT_V)) ASYNC(1)

    ISEND_WSET_OFFSET_V = ISEND_WSET_OFFSET(ISETW)
    ISEND_WSET_SIZE_V = ISEND_WSET_SIZE(ISETW)
    IF(PRESENT(PGP)) THEN
      !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(JK,JBLK,IFLD,JI) ASYNC(1)
      DO JFLD=1,ISEND_FIELD_COUNT_V
        DO JL=1,ISEND_WSET_SIZE_V
          JK = MOD(ISEND_WSET_OFFSET_V+JL-1,NPROMA)+1
          JBLK = (ISEND_WSET_OFFSET_V+JL-1)/NPROMA+1
          IFLD = IFLDA(JFLD)
          JI = (JFLD-1)*ISEND_WSET_SIZE_V+JL
          ZCOMBUFS(ICOMBUFS_OFFSET_V+JI) = PGP(JK,IFLD,JBLK)
        ENDDO
      ENDDO
    ELSE
      !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(JK,JBLK,IFLD,JI,IOFF) ASYNC(1)
      DO JFLD=1,ISEND_FIELD_COUNT_V
        DO JL=1,ISEND_WSET_SIZE_V
          JK = MOD(ISEND_WSET_OFFSET_V+JL-1,NPROMA)+1
          JBLK = (ISEND_WSET_OFFSET_V+JL-1)/NPROMA+1
          IFLD = IFLDA(JFLD)
          JI = ICOMBUFS_OFFSET_V+(JFLD-1)*ISEND_WSET_SIZE_V+JL
          IF(IFLD < PGP_INDICES(PGP_INDICES_UV+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_UV)
            PBOUND=UBOUND(PGPUV,2)
            ! TODO we could certainly reshape PGPXX arrays and we would simplify this
            ZCOMBUFS(JI) = PGPUV(JK,MOD(IOFF,PBOUND)+1,IOFF/PBOUND+1,JBLK)
          ELSEIF(IFLD < PGP_INDICES(PGP_INDICES_GP2+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_GP2)
            ZCOMBUFS(JI)  = PGP2(JK,IOFF+1,JBLK)
          ELSEIF(IFLD < PGP_INDICES(PGP_INDICES_GP3A+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_GP3A)
            PBOUND=UBOUND(PGP3A,2)
            ZCOMBUFS(JI) = PGP3A(JK,MOD(IOFF,PBOUND)+1,IOFF/PBOUND+1,JBLK)
          ELSEIF(IFLD < PGP_INDICES(PGP_INDICES_GP3B+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_GP3B)
            PBOUND=UBOUND(PGP3B,2)
            ZCOMBUFS(JI)= PGP3B(JK,MOD(IOFF,PBOUND)+1,IOFF/PBOUND+1,JBLK)
          ENDIF
       ENDDO
      ENDDO
    ENDIF
    !$ACC END DATA
  ENDDO

  IF (IRECV_COUNTS > 0) ALLOCATE(ZCOMBUFR(ICOMBUFR_OFFSET(IRECV_COUNTS+1)))
  !$ACC DATA IF(IRECV_COUNTS > 0) CREATE(ZCOMBUFR) ASYNC(1)

  !$ACC WAIT(1)

  CALL GSTATS(1602,1)

  IF (LSYNC_TRANS) THEN
    CALL GSTATS(430,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(430,1)
  ENDIF
  CALL GSTATS(411,0)

  IR=0

  !$ACC HOST_DATA USE_DEVICE(ZCOMBUFR,ZCOMBUFS)
  !  Receive loop.........................................................
  DO INR=1,IRECV_COUNTS
    IR=IR+1
    IPROC=IRECV_TO_PROC(INR)
    CALL MPI_IRECV(ZCOMBUFR(ICOMBUFR_OFFSET(INR)+1:ICOMBUFR_OFFSET(INR+1)),IRECVTOT(IPROC), &
      & TRGTOL_DTYPE,NPRCIDS(IPROC)-1,MTAGLG,MPL_COMM_OML(OML_MY_THREAD()),IREQ(IR),IERROR)
  ENDDO

  !....Send loop.........................................................
  DO INS=1,ISEND_COUNTS
    IR=IR+1
    ISEND=ISEND_TO_PROC(INS)
    CALL MPI_ISEND(ZCOMBUFS(ICOMBUFS_OFFSET(INS)+1:ICOMBUFS_OFFSET(INS+1)),ISENDTOT(ISEND), &
     & TRGTOL_DTYPE,NPRCIDS(ISEND)-1,MTAGLG,MPL_COMM_OML(OML_MY_THREAD()),IREQ(IR),IERROR)
  ENDDO
  !$ACC END HOST_DATA

  ! Copy local contribution
  IF(ISENDTOT(MYPROC) > 0 )THEN
    ! I have to send something to myself...

    ! Input is KF_GP fields. We find the resulting KF_FS fields.
    IFLDS = 0
    DO JFLD=1,KF_GP
      IF(KVSET(JFLD) == MYSETV .OR. KVSET(JFLD) == -1) THEN
        IFLDS = IFLDS+1
        IF(PRESENT(KPTRGP)) THEN
          IFLDA(IFLDS) = KPTRGP(JFLD)
        ELSE
          IFLDA(IFLDS) = JFLD
        ENDIF
      ENDIF
    ENDDO

    !$ACC DATA COPYIN(IFLDA(1:IFLDS)) ASYNC(1)

    ISEND_WSET_OFFSET_V = ISEND_WSET_OFFSET(MYSETW)
    ISEND_WSET_SIZE_V = ISEND_WSET_SIZE(MYSETW)
    IRECV_BUFR_TO_OUT_V = IRECV_BUFR_TO_OUT_OFFSET(MYPROC)
    CALL GSTATS(1601,0)
    IF(PRESENT(PGP)) THEN
      !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(JK,JBLK,IFLD,IPOS,IOFF) ASYNC(1)
      DO JFLD=1,KF_FS
        DO JL=1,ISEND_WSET_SIZE_V
          JK = MOD(ISEND_WSET_OFFSET_V+JL-1,NPROMA)+1
          JBLK = (ISEND_WSET_OFFSET_V+JL-1)/NPROMA+1
          IFLD = IFLDA(JFLD)
          IPOS = IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_V+JL,1)+ &
              & (JFLD-1)*IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_V+JL,2)+1
          PREEL_REAL(IPOS) = PGP(JK,IFLD,JBLK)
        ENDDO
      ENDDO
    ELSE
      !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(JK,JBLK,IFLD,IPOS,IOFF) ASYNC(1)
      DO JFLD=1,KF_FS
        DO JL=1,ISEND_WSET_SIZE_V
          JK = MOD(ISEND_WSET_OFFSET_V+JL-1,NPROMA)+1
          JBLK = (ISEND_WSET_OFFSET_V+JL-1)/NPROMA+1
          IFLD = IFLDA(JFLD)
          IPOS = IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_V+JL,1)+ &
              & (JFLD-1)*IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_V+JL,2)+1
          IF(IFLD < PGP_INDICES(PGP_INDICES_UV+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_UV)
            PBOUND=UBOUND(PGPUV,2)
            PREEL_REAL(IPOS) = PGPUV(JK,MOD(IOFF,PBOUND)+1,IOFF/PBOUND+1,JBLK)
          ELSEIF(IFLD < PGP_INDICES(PGP_INDICES_GP2+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_GP2)
            PREEL_REAL(IPOS) = PGP2(JK,IOFF+1,JBLK)
          ELSEIF(IFLD < PGP_INDICES(PGP_INDICES_GP3A+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_GP3A)
            PBOUND=UBOUND(PGP3A,2)
            PREEL_REAL(IPOS) = PGP3A(JK,MOD(IOFF,PBOUND)+1,IOFF/PBOUND+1,JBLK)
          ELSEIF(IFLD < PGP_INDICES(PGP_INDICES_GP3B+1)) THEN
            IOFF=IFLD-PGP_INDICES(PGP_INDICES_GP3B)
            PBOUND=UBOUND(PGP3B,2)
            PREEL_REAL(IPOS) = PGP3B(JK,MOD(IOFF,PBOUND)+1,IOFF/PBOUND+1,JBLK)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    CALL GSTATS(1601,1)

    !$ACC END DATA

  ENDIF


  IF(IR > 0) THEN
    CALL MPL_WAIT(KREQUEST=IREQ(1:IR), &
      & CDSTRING='TRGTOL_CUDAAWARE: WAIT FOR SENDS AND RECEIVES')
  ENDIF
  IF (LSYNC_TRANS) THEN
    CALL GSTATS(431,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(431,1)
  ENDIF
  CALL GSTATS(411,1)

  !$ACC EXIT DATA IF(ISEND_COUNTS > 0) DELETE(ZCOMBUFS)
  IF (ISEND_COUNTS > 0) DEALLOCATE(ZCOMBUFS)

  !  Unpack loop.........................................................

  CALL GSTATS(1603,0)


  DO INR=1,IRECV_COUNTS
    IPROC=IRECV_TO_PROC(INR)
    ILEN = IRECVTOT(IPROC)/KF_FS
    IRECV_BUFR_TO_OUT_V = IRECV_BUFR_TO_OUT_OFFSET(IPROC)
    ICOMBUFR_OFFSET_V = ICOMBUFR_OFFSET(INR)
    !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(II) ASYNC(1)
    DO JFLD=1,KF_FS
      DO JL=1,ILEN
        IPOS = IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_V+JL,1)+ &
            & (JFLD-1)*IRECV_BUFR_TO_OUT(IRECV_BUFR_TO_OUT_V+JL,2)+1
        PREEL_REAL(IPOS) = ZCOMBUFR(ICOMBUFR_OFFSET_V+JL+(JFLD-1)*ILEN)
      ENDDO
    ENDDO
  ENDDO

  !$ACC WAIT(1)
  CALL GSTATS(1603,1)

  !$ACC END DATA ! ZCOMBUFR
  !$ACC END DATA ! IRECV_BUFR_TO_OUT,PGPINDICES
  !$ACC END DATA !ZCOMBUFS (present)
  !$ACC END DATA !PGP3B
  !$ACC END DATA !PGP3A
  !$ACC END DATA !PGP2
  !$ACC END DATA !PGPUV
  !$ACC END DATA !PGP

  IF (IRECV_COUNTS > 0) DEALLOCATE(ZCOMBUFR)

  IF (LHOOK) CALL DR_HOOK('TRGTOL_CUDAAWARE',1,ZHOOK_HANDLE)

  END SUBROUTINE TRGTOL_CUDAAWARE

  END MODULE TRGTOL_MOD
