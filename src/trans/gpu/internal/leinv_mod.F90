#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEINV_MOD
  USE PARKIND_ECTRANS  ,ONLY : JPIM
  USE TPM_TRANS, ONLY: LEINV_CONFIG
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LEINV_STRIDES, LEINV
  PUBLIC :: SETUP_LEINV

  INTEGER(KIND=JPIM) :: A = 8 !Alignment

CONTAINS
  SUBROUTINE SETUP_LEINV(CONFIG)
    USE TPM_DISTR, ONLY: D
    USE TPM_DIM, ONLY: R
    USE TPM_GEOMETRY, ONLY: G
    USE PARKIND_ECTRANS, ONLY: JPIM

    IMPLICIT NONE

    TYPE(LEINV_CONFIG), INTENT(OUT) :: CONFIG

    INTEGER(KIND=JPIM) :: KM, KMLOC, N_OFFSET, K_OFFSET

    ALLOCATE(CONFIG%OFFSETS_N(D%NUMP+1))
    ALLOCATE(CONFIG%OFFSETS_K(D%NUMP+1))

    N_OFFSET = 0
    K_OFFSET = 0
    DO KMLOC=1,D%NUMP
      KM = D%MYMS(KMLOC)
      CONFIG%OFFSETS_K(KMLOC) = K_OFFSET
      CONFIG%OFFSETS_N(KMLOC) = N_OFFSET

      !KM=0 is transformed in double precision, no need to store here
      IF (KM /= 0) THEN
        ! K_OFFSET takes the max of the two GEMMs
        K_OFFSET = K_OFFSET + ALIGN((R%NSMAX-KM+3)/2,A)
        N_OFFSET = N_OFFSET + ALIGN(G%NDGLU(KM),A)
      ENDIF
    ENDDO
    CONFIG%OFFSETS_K(D%NUMP+1) = K_OFFSET
    CONFIG%OFFSETS_N(D%NUMP+1) = N_OFFSET

    !$ACC ENTER DATA COPYIN(CONFIG,CONFIG%OFFSETS_K,CONFIG%OFFSETS_N)
  END SUBROUTINE
  SUBROUTINE LEINV_STRIDES(KF_LEG,IOUT_STRIDES0,IOUT_SIZE,IIN_STRIDES0,IIN_SIZE,&
                           IOUT0_STRIDES0,IOUT0_SIZE,IIN0_STRIDES0,IIN0_SIZE)
    USE PARKIND_ECTRANS  ,ONLY : JPIM, JPRBT, JPRD
    USE TPM_DIM         ,ONLY : R
    USE TPM_DISTR, ONLY: D
    USE TPM_TRANS, ONLY: LEINV_CONF

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN)  :: KF_LEG

    INTEGER(KIND=JPIM), OPTIONAL :: IOUT_STRIDES0, IOUT_SIZE
    INTEGER(KIND=JPIM), OPTIONAL :: IIN_STRIDES0, IIN_SIZE
    INTEGER(KIND=JPIM), OPTIONAL :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM), OPTIONAL :: IIN0_STRIDES0, IIN0_SIZE


    IF (PRESENT(IOUT_STRIDES0)) &
      IOUT0_STRIDES0 = ALIGN(KF_LEG,A)
    IF (PRESENT(IOUT0_SIZE)) &
      IOUT0_SIZE = IOUT0_STRIDES0 * ALIGN(R%NDGNH,A)
    IF (PRESENT(IIN_STRIDES0)) &
      IIN_STRIDES0 = ALIGN(2*KF_LEG,A)
    IF (PRESENT(IIN_SIZE)) &
      IIN_SIZE = IIN_STRIDES0*LEINV_CONF%OFFSETS_K(D%NUMP+1)
    IF (PRESENT(IOUT0_STRIDES0)) &
      IOUT_STRIDES0 = ALIGN(2*KF_LEG,A)
    IF (PRESENT(IOUT_SIZE)) &
      IOUT_SIZE = IOUT_STRIDES0*LEINV_CONF%OFFSETS_N(D%NUMP+1)
    IF (PRESENT(IIN0_STRIDES0)) &
      IIN0_STRIDES0 = ALIGN(KF_LEG,A)
    IF (PRESENT(IIN0_SIZE)) &
      IIN0_SIZE = IIN0_STRIDES0 * ALIGN(MAX((R%NTMAX+2)/2,(R%NTMAX+3)/2),A)
  END SUBROUTINE

  SUBROUTINE LEINV(PIA,ZINP,ZINP0,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,KF_LEG)
    !**** *LEINV* - Inverse Legendre transform.

    !     Purpose.
    !     --------
    !        Inverse Legendre tranform of all variables(kernel).

    !**   Interface.
    !     ----------
    !        CALL LEINV(...)

    !        Explicit arguments :  KM - zonal wavenumber (input-c)
    !        --------------------  KFC - number of fields to tranform (input-c)
    !                              PIA - spectral fields
    !                              for zonal wavenumber KM (input)

    !        Implicit arguments :  None.
    !        --------------------

    !     Method.
    !     -------

    !     Externals.
    !     ----------

    !     Reference.
    !     ----------
    !        ECMWF Research Department documentation of the IFS

    !     Author.
    !     -------
    !      Nils Wedi + Mats Hamrud + George Modzynski
    !
    !     Modifications.
    !     --------------
    !        J.Hague : Oct 2012 DR_HOOK round calls to DGEMM:
    !      F. Vana  05-Mar-2015  Support for single precision
    !     ------------------------------------------------------------------

    USE TPM_GEN         ,ONLY : LSYNC_TRANS
    USE PARKIND_ECTRANS  ,ONLY : JPIM,JPRB,  JPRBT, JPRD
    USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
    USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX, R_NDGL
    USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
    USE TPM_FIELDS      ,ONLY : ZAA,ZAS,ZAA0,ZAS0,KMLOC0
    USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,MYPROC
    USE CUDA_GEMM_BATCHED_MOD
    USE MPL_MODULE      ,ONLY : MPL_BARRIER
    USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX
    USE TPM_TRANS, ONLY: LEINV_CONF

    IMPLICIT NONE

    REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN)  :: KF_LEG
    REAL(KIND=JPRBT), INTENT(OUT) :: ZINP(:), ZOUTS(:), ZOUTA(:)
    REAL(KIND=JPRD), INTENT(OUT) :: ZINP0(:), ZOUTS0(:), ZOUTA0(:)

    !     LOCAL
    INTEGER(KIND=JPIM)  :: KS(D_NUMP), NS(D_NUMP), AOFFSETS(D_NUMP), BOFFSETS(D_NUMP), COFFSETS(D_NUMP)
    INTEGER(KIND=JPIM)  :: KM, KMLOC, IA, IS, ISL, J1, JGL, JK, J
    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0, IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IIN_STRIDES0, IIN_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM)  :: IIN0_STRIDES0, IIN0_SIZE

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


    !*       1.1      PREPARATIONS.
    IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)

    !     ------------------------------------------------------------------

    !*       1.       PERFORM LEGENDRE TRANFORM.
    !                 --------------------------

    !*       1.1      PREPARATIONS.

    CALL LEINV_STRIDES(KF_LEG,IOUT_STRIDES0,IOUT_SIZE,IIN_STRIDES0,IIN_SIZE,&
                       IOUT0_STRIDES0,IOUT0_SIZE,IIN0_STRIDES0,IIN0_SIZE)


    !$ACC DATA PRESENT(D,D_MYMS,G,G_NDGLU,R) &
    !$ACC&     PRESENT(ZINP,ZOUTS,ZOUTA,ZINP0,ZOUTS0,ZOUTA0) &
    !$ACC&     PRESENT(ZAA,ZAS,PIA) &
    !$ACC&     PRESENT(D_MYMS,G_NDGLU)

    IF (KMLOC0 > 0) THEN
      print*,'computing m=0 in double precision'
    ENDIF

    ! READ 2:NSMAX+3

    !IF KM=0 and NSMAX is 6:
    !    IA=1
    !    DO=1,6/2+1 ... 1..4
    !       PIA_2=1+1+(J-1)*2 ...2+(0..3)*2 .... 2,4,6,8
    !IF KM=0 and NSMAX is 7:
    !    IA=2
    !    DO=1,7/2+1 ... 1..4
    !       PIA_2=2+1+(1..4-1)*2 ...3+(0..3)*2 .... 3,5,7,9

    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IA,J) DEFAULT(NONE) ASYNC(1)
    DO KMLOC=1,D_NUMP
      DO JK=1,2*KF_LEG
        KM =  D_MYMS(KMLOC)
        IA  = 1+MOD(R_NSMAX-KM+2,2)
        IF(KM /= 0)THEN
          !$ACC LOOP SEQ
          DO J=1,(R_NSMAX-KM+2)/2
            ZINP(JK+(J-1)*IIN_STRIDES0+LEINV_CONF%OFFSETS_K(KMLOC)*IIN_STRIDES0)=PIA(JK,IA+1+(J-1)*2,KMLOC)
          ENDDO
        ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
          ! every other field is sufficient because Im(KM=0) == 0
          !$ACC LOOP SEQ
          DO J=1,(R_NSMAX+2)/2
            ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0) = PIA(JK,IA+1+(J-1)*2,KMLOC)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    IF (LSYNC_TRANS) THEN
      !$ACC WAIT(1)
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(424,0)

    IF (KMLOC0 > 0) THEN
      ! compute m=0 in double precision:
      CALL CUDA_GEMM_BATCHED( &
        & CUBLAS_OP_N, CUBLAS_OP_T, &
        & KF_LEG, G%NDGLU(0), (R%NSMAX+2)/2, &
        & 1.0_JPRD, &
        & ZINP0, IIN0_STRIDES0, 0, &
        & ZAA0, SIZE(ZAA0,1), 0, &
        & 0.0_JPRD, &
        & ZOUTA0, IOUT0_STRIDES0, 0, &
        & 1, STREAM=1_C_LONG)
    ENDIF

    DO KMLOC=1,D_NUMP
      KM = D_MYMS(KMLOC)
      KS(KMLOC) = (R%NSMAX-KM+2)/2
      NS(KMLOC) = G%NDGLU(KM)
      AOFFSETS(KMLOC) = IIN_STRIDES0*LEINV_CONF%OFFSETS_K(KMLOC)
      BOFFSETS(KMLOC) = SIZE(ZAA,1)*SIZE(ZAA,2)*(KMLOC-1)
      COFFSETS(KMLOC) = IOUT_STRIDES0*LEINV_CONF%OFFSETS_N(KMLOC)
    ENDDO
    IF(KMLOC0 > 0) THEN
      NS(KMLOC0) = 0
      KS(KMLOC0) = 0
    ENDIF
    CALL CUDA_GEMM_BATCHED( &
        & 11, & ! unique identifier
        & CUBLAS_OP_N, CUBLAS_OP_T, &
        & 2*KF_LEG, NS(:), KS(:), &
        & 1.0_JPRBT, &
        & ZINP, IIN_STRIDES0, AOFFSETS, &
        & ZAA, SIZE(ZAA,1), BOFFSETS, &
        & 0.0_JPRBT, &
        & ZOUTA, IOUT_STRIDES0, COFFSETS, &
        & D_NUMP, STREAM=1_C_LONG)

    IF (LSYNC_TRANS) THEN
      !$ACC WAIT(1)
      CALL GSTATS(444,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(444,1)
    ENDIF
    CALL GSTATS(424,1)

    ! 2. +++++++++++++ symmetric
    !IF KM=0 and NSMAX is 6:
    !    IS=2
    !    DO=1,4
    !       PIA_2=2+1+(0..3)*2 ... 3+(0..3)*2 ... 3,5,7,9
    !IF KM=0 and NSMAX is 7:
    !    IS=1
    !    DO=1,5
    !       PIA_2=1+1+(1..5-1)*2 ...2+(0..4)*2 .... 2,4,6,8,10

    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IS,J) DEFAULT(NONE) ASYNC(1)
    DO KMLOC=1,D_NUMP
      DO JK=1,2*KF_LEG
        KM =  D_MYMS(KMLOC)
        IS  = 1+MOD(R_NSMAX-KM+1,2)
        IF(KM /= 0) THEN
          !$ACC LOOP SEQ
          DO J=1,(R_NSMAX-KM+3)/2
            ZINP(JK+(J-1)*IIN_STRIDES0+LEINV_CONF%OFFSETS_K(KMLOC)*IIN_STRIDES0)=PIA(JK,IS+1+(J-1)*2,KMLOC)
          ENDDO
        ELSEIF (MOD((JK-1),2) == 0) THEN
          !$ACC LOOP SEQ
          DO J=1,(R_NSMAX+3)/2
            ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0) = PIA(JK,IS+1+(J-1)*2,KMLOC)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    IF (LSYNC_TRANS) THEN
      !$ACC WAIT(1)
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(424,0)

    IF (KMLOC0 > 0) THEN
      CALL CUDA_GEMM_BATCHED( &
        & CUBLAS_OP_N, CUBLAS_OP_T, &
        & KF_LEG, G%NDGLU(0), (R%NSMAX+3)/2, &
        & 1.0_JPRD, &
        & ZINP0, IIN0_STRIDES0, 0, &
        & ZAS0, SIZE(ZAS0,1), 0, &
        & 0.0_JPRD, &
        & ZOUTS0, IOUT0_STRIDES0, 0, &
        & 1, STREAM=1_C_LONG)
    ENDIF

    DO KMLOC=1,D_NUMP
      KM = D_MYMS(KMLOC)
      KS(KMLOC) = (R%NSMAX-KM+3)/2
      NS(KMLOC) = G%NDGLU(KM)
      AOFFSETS(KMLOC) = IIN_STRIDES0*LEINV_CONF%OFFSETS_K(KMLOC)
      BOFFSETS(KMLOC) = SIZE(ZAS,1)*SIZE(ZAS,2)*(KMLOC-1)
      COFFSETS(KMLOC) = IOUT_STRIDES0*LEINV_CONF%OFFSETS_N(KMLOC)
    ENDDO
    IF(KMLOC0 > 0) THEN
      NS(KMLOC0) = 0
      KS(KMLOC0) = 0
    ENDIF
    CALL CUDA_GEMM_BATCHED( &
      & 12, & ! unique identifier
      & CUBLAS_OP_N, CUBLAS_OP_T, &
      & 2*KF_LEG, NS(:), KS(:), &
      & 1.0_JPRBT, &
      & ZINP, IIN_STRIDES0, AOFFSETS, &
      & ZAS, SIZE(ZAS,1), BOFFSETS, &
      & 0.0_JPRBT, &
      & ZOUTS, IOUT_STRIDES0, COFFSETS, &
      & D_NUMP, STREAM=1_C_LONG)
    IF (LSYNC_TRANS) THEN
      !$ACC WAIT(1)
      CALL GSTATS(444,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(444,1)
    ENDIF
    CALL GSTATS(424,1)

    !$ACC WAIT(1)

    !$ACC END DATA

    IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
    !     ------------------------------------------------------------------
  END SUBROUTINE LEINV
END MODULE LEINV_MOD
