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

MODULE LEINV_MOD
  USE PARKIND_ECTRANS  ,ONLY : JPIM
  USE BUFFERED_ALLOCATOR_MOD
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LEINV_STRIDES, LEINV

  INTEGER(KIND=JPIM) :: A = 8 !Alignment

CONTAINS
  SUBROUTINE LEINV_STRIDES(KF_LEG,IOUT_STRIDES0,IOUT_SIZE,IIN_STRIDES0,IIN_SIZE,&
                           IOUT0_STRIDES0,IOUT0_SIZE,IIN0_STRIDES0,IIN0_SIZE)
    USE PARKIND_ECTRANS  ,ONLY : JPIM, JPRBT, JPRD
    USE TPM_DIM         ,ONLY : R
    USE TPM_DISTR, ONLY: D,D_OFFSETS_GEMM1,D_OFFSETS_GEMM2

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
      IIN_SIZE = IIN_STRIDES0*D_OFFSETS_GEMM2(D%NUMP+1)
    IF (PRESENT(IOUT0_STRIDES0)) &
      IOUT_STRIDES0 = ALIGN(2*KF_LEG,A)
    IF (PRESENT(IOUT_SIZE)) &
      IOUT_SIZE = IOUT_STRIDES0*D_OFFSETS_GEMM1(D%NUMP+1)
    IF (PRESENT(IIN0_STRIDES0)) &
      IIN0_STRIDES0 = ALIGN(KF_LEG,A)
    IF (PRESENT(IIN0_SIZE)) &
      IIN0_SIZE = IIN0_STRIDES0 * ALIGN(MAX((R%NTMAX+2)/2,(R%NTMAX+3)/2),A)
  END SUBROUTINE LEINV_STRIDES

  SUBROUTINE LEINV(ALLOCATOR,PIA,ZINP,ZINP0,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,KF_LEG)
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
    USE TPM_DIM         ,ONLY : R_NDGNH,R_NSMAX, R_NDGL
    USE TPM_GEOMETRY    ,ONLY : G_NDGLU
    USE TPM_FIELDS      ,ONLY : ZAA,ZAS,ZAA0,ZAS0,KMLOC0
    USE TPM_DISTR       ,ONLY : D_NUMP,D_MYMS,MYPROC,D_OFFSETS_GEMM1,D_OFFSETS_GEMM2
    USE HICBLAS_MOD     ,ONLY : HIP_GEMM_BATCHED, HIP_DGEMM_BATCHED_OVERLOAD, &
 &                              HIP_DGEMM_GROUPED_OVERLOAD, HIP_SGEMM_GROUPED_OVERLOAD
#ifdef TRANS_SINGLE
#define HIP_GEMM HIP_SGEMM_GROUPED_OVERLOAD
#else
#define HIP_GEMM HIP_DGEMM_GROUPED_OVERLOAD
#endif

    USE, INTRINSIC :: ISO_C_BINDING
    USE MPL_MODULE      ,ONLY : MPL_BARRIER
    USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

    IMPLICIT NONE

    REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN)  :: KF_LEG
    REAL(KIND=JPRBT), INTENT(OUT) :: ZINP(:), ZOUTS(:), ZOUTA(:)
    REAL(KIND=JPRD), INTENT(OUT) :: ZINP0(:), ZOUTS0(:), ZOUTA0(:)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR

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


#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC DATA PRESENT(D_MYMS,D_NUMP,G_NDGLU) &
    !$ACC&     PRESENT(ZINP,ZOUTS,ZOUTA,ZINP0,ZOUTS0,ZOUTA0) &
    !$ACC&     PRESENT(ZAA,ZAS,PIA) &
    !$ACC&     PRESENT(R_NSMAX,G_NDGLU,D_OFFSETS_GEMM2)
#endif

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

#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IA,J) FIRSTPRIVATE(KF_LEG,IIN_STRIDES0,IIN0_STRIDES0) DEFAULT(NONE) ASYNC(1)
#endif
    DO KMLOC=1,D_NUMP
      DO JK=1,2*KF_LEG
        KM =  D_MYMS(KMLOC)
        IA  = 1+MOD(R_NSMAX-KM+2,2)
        IF(KM /= 0)THEN
#ifdef OMPGPU
#endif
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX-KM+2)/2
            ZINP(JK+(J-1)*IIN_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IIN_STRIDES0)=PIA(JK,IA+1+(J-1)*2,KMLOC)
          ENDDO
        ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
          ! every other field is sufficient because Im(KM=0) == 0
#ifdef OMPGPU
#endif
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX+2)/2
            ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0) = PIA(JK,IA+1+(J-1)*2,KMLOC)
          ENDDO
        ENDIF
      ENDDO
    ENDDO


    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(424,0)

    IF (KMLOC0 > 0) THEN
      ! compute m=0 in double precision
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAA0,ZINP0,ZOUTA0)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAA0,ZINP0,ZOUTA0)
#endif
      CALL HIP_DGEMM_BATCHED_OVERLOAD( &
        & 'N', 'T', &
        & KF_LEG, G_NDGLU(0), (R_NSMAX+2)/2, &
        & 1.0_JPRD, &
        & ZINP0, IIN0_STRIDES0, 0, &
        & ZAA0, SIZE(ZAA0,1), 0, &
        & 0.0_JPRD, &
        & ZOUTA0, IOUT0_STRIDES0, 0, &
        & 1, STREAM=1_C_INT, ALLOC=ALLOCATOR%PTR)
#ifdef OMPGPU
      !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
      !$ACC END HOST_DATA
#endif
   ENDIF



    DO KMLOC=1,D_NUMP
      KM = D_MYMS(KMLOC)
      KS(KMLOC) = (R_NSMAX-KM+2)/2
      NS(KMLOC) = G_NDGLU(KM)
      AOFFSETS(KMLOC) = IIN_STRIDES0*D_OFFSETS_GEMM2(KMLOC)
      BOFFSETS(KMLOC) = SIZE(ZAA,1)*SIZE(ZAA,2)*(KMLOC-1)
      COFFSETS(KMLOC) = IOUT_STRIDES0*D_OFFSETS_GEMM1(KMLOC)
    ENDDO
    IF(KMLOC0 > 0) THEN
      NS(KMLOC0) = 0
      KS(KMLOC0) = 0
    ENDIF
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAA,ZINP,ZOUTA)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAA,ZINP,ZOUTA)
#endif
    CALL HIP_GEMM( &
        & 11, & ! unique identifier
        & 'N', 'T', &
        & 2*KF_LEG, NS(:), KS(:), &
        & 1.0_JPRBT, &
        & ZINP, IIN_STRIDES0, AOFFSETS, &
        & ZAA, SIZE(ZAA,1), BOFFSETS, &
        & 0.0_JPRBT, &
        & ZOUTA, IOUT_STRIDES0, COFFSETS, &
        & D_NUMP, STREAM=1_C_INT, ALLOC=ALLOCATOR%PTR)
#ifdef OMPGPU
      !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
      !$ACC END HOST_DATA
#endif

    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
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

#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IS,J) FIRSTPRIVATE(KF_LEG,IIN_STRIDES0,IIN0_STRIDES0) DEFAULT(NONE) ASYNC(1)
#endif
    DO KMLOC=1,D_NUMP
      DO JK=1,2*KF_LEG
        KM =  D_MYMS(KMLOC)
        IS  = 1+MOD(R_NSMAX-KM+1,2)
        IF(KM /= 0) THEN
#ifdef OMPGPU
#endif
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX-KM+3)/2
            ZINP(JK+(J-1)*IIN_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IIN_STRIDES0)=PIA(JK,IS+1+(J-1)*2,KMLOC)
          ENDDO
        ELSEIF (MOD((JK-1),2) == 0) THEN
#ifdef OMPGPU
#endif
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX+3)/2
            ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0) = PIA(JK,IS+1+(J-1)*2,KMLOC)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(424,0)

    IF (KMLOC0 > 0) THEN
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAS0,ZINP0,ZOUTS0)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAS0,ZINP0,ZOUTS0)
#endif
      CALL HIP_DGEMM_BATCHED_OVERLOAD( &
        & 'N', 'T', &
        & KF_LEG, G_NDGLU(0), (R_NSMAX+3)/2, &
        & 1.0_JPRD, &
        & ZINP0, IIN0_STRIDES0, 0, &
        & ZAS0, SIZE(ZAS0,1), 0, &
        & 0.0_JPRD, &
        & ZOUTS0, IOUT0_STRIDES0, 0, &
        & 1, STREAM=1_C_INT, ALLOC=ALLOCATOR%PTR)
#ifdef OMPGPU
      !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
      !$ACC END HOST_DATA
#endif
    ENDIF

    DO KMLOC=1,D_NUMP
      KM = D_MYMS(KMLOC)
      KS(KMLOC) = (R_NSMAX-KM+3)/2
      NS(KMLOC) = G_NDGLU(KM)
      AOFFSETS(KMLOC) = IIN_STRIDES0*D_OFFSETS_GEMM2(KMLOC)
      BOFFSETS(KMLOC) = SIZE(ZAS,1)*SIZE(ZAS,2)*(KMLOC-1)
      COFFSETS(KMLOC) = IOUT_STRIDES0*D_OFFSETS_GEMM1(KMLOC)
    ENDDO
    IF(KMLOC0 > 0) THEN
      NS(KMLOC0) = 0
      KS(KMLOC0) = 0
    ENDIF
#ifdef OMPGPU
    !$OMP TARGET DATA USE_DEVICE_PTR(ZAS,ZINP,ZOUTS)
#endif
#ifdef ACCGPU
    !$ACC HOST_DATA USE_DEVICE(ZAS,ZINP,ZOUTS)
#endif
    CALL HIP_GEMM( &
      & 12, & ! unique identifier
      & 'N', 'T', &
      & 2*KF_LEG, NS(:), KS(:), &
      & 1.0_JPRBT, &
      & ZINP, IIN_STRIDES0, AOFFSETS, &
      & ZAS, SIZE(ZAS,1), BOFFSETS, &
      & 0.0_JPRBT, &
      & ZOUTS, IOUT_STRIDES0, COFFSETS, &
      & D_NUMP, STREAM=1_C_INT, ALLOC=ALLOCATOR%PTR)
#ifdef OMPGPU
    !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
    !$ACC END HOST_DATA
#endif
    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(444,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(444,1)
    ENDIF
    CALL GSTATS(424,1)

#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC WAIT(1)

    !$ACC END DATA
#endif

    IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
    !     ------------------------------------------------------------------
  END SUBROUTINE LEINV
END MODULE LEINV_MOD
