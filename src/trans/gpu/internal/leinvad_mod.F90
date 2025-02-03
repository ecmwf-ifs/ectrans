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

MODULE LEINVAD_MOD
  USE PARKIND_ECTRANS,        ONLY: JPIM, JPRB, JPRBT, JPRD, JPIB
  USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR
  USE LEINV_MOD,              ONLY: LEINV_STRIDES
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LEINVAD

  INTEGER(KIND=JPIM) :: A = 8 !Alignment

CONTAINS

  SUBROUTINE LEINVAD(ALLOCATOR,PIA,ZINP,ZINP0,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,KF_LEG)
    !**** *LEINVAD* - Adjoint of inverse Legendre transform.

    !     Purpose.
    !     --------
    !        Adjoint of inverse Legendre tranform of all variables(kernel).

    !**   Interface.
    !     ----------
    !        CALL LEINVAD(...)

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

    USE TPM_GEN,                     ONLY: LSYNC_TRANS, NOUT, NCUR_RESOL
    USE YOMHOOK,                     ONLY: LHOOK, DR_HOOK, JPHOOK
    USE TPM_DIM,                     ONLY: R
    USE TPM_GEOMETRY,                ONLY: G
    USE TPM_FIELDS_GPU,              ONLY: FG
    USE TPM_DISTR,                   ONLY: D
    USE HICBLAS_MOD,                 ONLY: HIP_DGEMM_BATCHED_OVERLOAD, &
      &                                    HIP_DGEMM_GROUPED_OVERLOAD, HIP_SGEMM_GROUPED_OVERLOAD
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
    USE MPL_MODULE,                  ONLY: MPL_BARRIER,MPL_ALL_MS_COMM
    USE TPM_STATS,                   ONLY: GSTATS => GSTATS_NVTX
#ifdef TRANS_SINGLE
#define HIP_GEMM HIP_SGEMM_GROUPED_OVERLOAD
#else
#define HIP_GEMM HIP_DGEMM_GROUPED_OVERLOAD
#endif

    IMPLICIT NONE

    REAL(KIND=JPRB),    INTENT(INOUT)  :: PIA(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN)  :: KF_LEG
    REAL(KIND=JPRBT), INTENT(IN) :: ZINP(:), ZOUTS(:), ZOUTA(:)
    REAL(KIND=JPRD), INTENT(IN) :: ZINP0(:), ZOUTS0(:), ZOUTA0(:)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR

    !     LOCAL
    INTEGER(KIND=JPIM)  :: KS(D%NUMP), NS(D%NUMP)
    INTEGER(KIND=JPIB)  :: AOFFSETS(D%NUMP), BOFFSETS(D%NUMP), COFFSETS(D%NUMP)
    INTEGER(KIND=JPIM)  :: KM, KMLOC, IA, IS, ISL, J1, JGL, JK, J, IMLOC0(1)
    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0
    INTEGER(KIND=JPIB)  :: IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IIN_STRIDES0
    INTEGER(KIND=JPIB)  :: IIN_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM)  :: IIN0_STRIDES0, IIN0_SIZE

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    ASSOCIATE(D_NUMP=>D%NUMP, R_NSMAX=>R%NSMAX, G_NDGLU=>G%NDGLU, D_MYMS=>D%MYMS, D_OFFSETS_GEMM1=>D%OFFSETS_GEMM1,&
        D_OFFSETS_GEMM2=>D%OFFSETS_GEMM2, &
        ZAA=>FG%ZAA, ZAS=>FG%ZAS, ZAA0=>FG%ZAA0, ZAS0=>FG%ZAS0)

    !*       1.1      PREPARATIONS.
    IF (LHOOK) CALL DR_HOOK('LEINVAD_DGEMM',0,ZHOOK_HANDLE)

    !     ------------------------------------------------------------------

    !*       1.       PERFORM LEGENDRE TRANFORM.
    !                 --------------------------

    !*       1.1      PREPARATIONS.

    CALL LEINV_STRIDES(KF_LEG,IOUT_STRIDES0,IOUT_SIZE,IIN_STRIDES0,IIN_SIZE,&
                       IOUT0_STRIDES0,IOUT0_SIZE,IIN0_STRIDES0,IIN0_SIZE)


#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC DATA PRESENT(D,D_MYMS,D_NUMP) &
    !$ACC&     PRESENT(ZINP,ZOUTS,ZOUTA,ZINP0,ZOUTS0,ZOUTA0) &
    !$ACC&     PRESENT(ZAA,ZAS,PIA) &
    !$ACC&     PRESENT(R,R_NSMAX,D_OFFSETS_GEMM2)
#endif

    ! READ 2:NSMAX+3

    ! 1. +++++++++++++ anti-symmetric
    !IF KM=0 and NSMAX is 6:
    !    IA=1
    !    DO=1,6/2+1 ... 1..4
    !       PIA_2=1+1+(J-1)*2 ...2+(0..3)*2 .... 2,4,6,8
    !IF KM=0 and NSMAX is 7:
    !    IA=2
    !    DO=1,7/2+1 ... 1..4
    !       PIA_2=2+1+(1..4-1)*2 ...3+(0..3)*2 .... 3,5,7,9

    CALL GSTATS(424,0)

    IMLOC0 = FINDLOC(D_MYMS,0)
    IF (IMLOC0(1) > 0) THEN
      ! compute m=0 in double precision
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAA0,ZINP0,ZOUTA0)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAA0,ZINP0,ZOUTA0)
#endif
      CALL HIP_DGEMM_BATCHED_OVERLOAD( &
        & 'N', 'N', &
        & KF_LEG, (R_NSMAX+2)/2, G_NDGLU(0), &
        & 1.0_JPRD, &
        & ZOUTA0, IOUT0_STRIDES0, 0, &
        & ZAA0, SIZE(ZAA0,1), 0, &
        & 0.0_JPRD, &
        & ZINP0, IIN0_STRIDES0, 0, &
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
      BOFFSETS(KMLOC) = D%OFFSETS_GEMM_MATRIX(KMLOC)
      COFFSETS(KMLOC) = IOUT_STRIDES0*D_OFFSETS_GEMM1(KMLOC)
    ENDDO
    IF(IMLOC0(1) > 0) THEN
      NS(IMLOC0(1)) = 0
      KS(IMLOC0(1)) = 0
    ENDIF
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAA,ZINP,ZOUTA)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAA,ZINP,ZOUTA)
#endif
    CALL HIP_GEMM( &
        & NCUR_RESOL, 11, & ! unique identifier
        & 'N', 'N', &
        & 2*KF_LEG, KS(:), NS(:), &
        & 1.0_JPRBT, &
        & ZOUTA, IOUT_STRIDES0, COFFSETS, &
        & ZAA, D%LEGENDRE_MATRIX_STRIDES, BOFFSETS, &
        & 0.0_JPRBT, &
        & ZINP, IIN_STRIDES0, AOFFSETS, &
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
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(444,1)
    ENDIF
    CALL GSTATS(424,1)

#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IA,J) &
    !$ACC& FIRSTPRIVATE(KF_LEG,IIN_STRIDES0,IIN0_STRIDES0) DEFAULT(NONE) &
#ifdef _CRAYFTN
    !$ACC&
#else
    !$ACC& ASYNC(1)
#endif
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
            PIA(JK,IA+1+(J-1)*2,KMLOC) = ZINP(JK+(J-1)*IIN_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IIN_STRIDES0)
          ENDDO
!           ! those are only needed with tensor cores (zinp might contain NaNs!)
! #if defined(USE_CUTLASS) && defined(USE_CUTLASS_3XTF32)
!           !$ACC LOOP SEQ
!           DO J=(R_NSMAX-KM+2)/2+1,ALIGN((R_NSMAX-KM+2)/2,A)
!             ZINP(JK+(J-1)*IIN_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IIN_STRIDES0)=0
!           ENDDO
! #endif
        ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
          ! every other field is sufficient because Im(KM=0) == 0
#ifdef OMPGPU
#endif
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX+2)/2
            PIA(JK,IA+1+(J-1)*2,KMLOC) = ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0)
          ENDDO
!           ! those are only needed with tensor cores (zinp might contain NaNs!)
! #if defined(USE_CUTLASS) && defined(USE_CUTLASS_3XTF32)
!           !$ACC LOOP SEQ
!           DO J=(R_NSMAX+2)/2+1,ALIGN((R_NSMAX+2)/2,A)
!             ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0) = 0
!           ENDDO
! #endif
        ENDIF
      ENDDO
    ENDDO


    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF


    ! 2. +++++++++++++ symmetric
    !IF KM=0 and NSMAX is 6:
    !    IS=2
    !    DO=1,4
    !       PIA_2=2+1+(0..3)*2 ... 3+(0..3)*2 ... 3,5,7,9
    !IF KM=0 and NSMAX is 7:
    !    IS=1
    !    DO=1,5
    !       PIA_2=1+1+(1..5-1)*2 ...2+(0..4)*2 .... 2,4,6,8,10

    CALL GSTATS(424,0)

    IF (IMLOC0(1) > 0) THEN
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAS0,ZINP0,ZOUTS0)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAS0,ZINP0,ZOUTS0)
#endif
      CALL HIP_DGEMM_BATCHED_OVERLOAD( &
        & 'N', 'N', &
        & KF_LEG, (R_NSMAX+3)/2, G_NDGLU(0), &
        & 1.0_JPRD, &
        & ZOUTS0, IOUT0_STRIDES0, 0, &
        & ZAS0, SIZE(ZAS0,1), 0, &
        & 0.0_JPRD, &
        & ZINP0, IIN0_STRIDES0, 0, &
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
      BOFFSETS(KMLOC) = D%OFFSETS_GEMM_MATRIX(KMLOC)
      COFFSETS(KMLOC) = IOUT_STRIDES0*D_OFFSETS_GEMM1(KMLOC)
    ENDDO
    IF(IMLOC0(1) > 0) THEN
      NS(IMLOC0(1)) = 0
      KS(IMLOC0(1)) = 0
    ENDIF
#ifdef OMPGPU
    !$OMP TARGET DATA USE_DEVICE_PTR(ZAS,ZINP,ZOUTS)
#endif
#ifdef ACCGPU
    !$ACC HOST_DATA USE_DEVICE(ZAS,ZINP,ZOUTS)
#endif
    CALL HIP_GEMM( &
      & NCUR_RESOL, 12, & ! unique identifier
      & 'N', 'N', &
      & 2*KF_LEG, KS(:), NS(:), &
      & 1.0_JPRBT, &
      & ZOUTS, IOUT_STRIDES0, COFFSETS, &
      & ZAS, D%LEGENDRE_MATRIX_STRIDES, BOFFSETS, &
      & 0.0_JPRBT, &
      & ZINP, IIN_STRIDES0, AOFFSETS, &
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
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(444,1)
    ENDIF
    CALL GSTATS(424,1)

#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IS,J) &
    !$ACC& FIRSTPRIVATE(KF_LEG,IIN_STRIDES0,IIN0_STRIDES0) DEFAULT(NONE) &
#ifndef _CRAYFTN
    !$ACC& ASYNC(1)
#else
    !$ACC&
#endif
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
            PIA(JK,IS+1+(J-1)*2,KMLOC) = ZINP(JK+(J-1)*IIN_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IIN_STRIDES0)
          ENDDO
!           ! those are only needed with tensor cores (zinp might contain NaNs!)
! #if defined(USE_CUTLASS) && defined(USE_CUTLASS_3XTF32)
!           !$ACC LOOP SEQ
!           DO J=(R_NSMAX-KM+3)/2+1,ALIGN((R_NSMAX-KM+3)/2,A)
!             ZINP(JK+(J-1)*IIN_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IIN_STRIDES0)=0
!           ENDDO
! #endif
        ELSEIF (MOD((JK-1),2) == 0) THEN
#ifdef OMPGPU
#endif
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX+3)/2
            PIA(JK,IS+1+(J-1)*2,KMLOC) = ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0)
          ENDDO
!           ! those are only needed with tensor cores (zinp might contain NaNs!)
! #if defined(USE_CUTLASS) && defined(USE_CUTLASS_3XTF32)
!           !$ACC LOOP SEQ
!           DO J=(R_NSMAX+3)/2+1,ALIGN((R_NSMAX+3)/2,A)
!             ZINP0((JK-1)/2+1+(J-1)*IIN0_STRIDES0) = 0
!           ENDDO
! #endif
        ENDIF
      ENDDO
    ENDDO

    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF


#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC WAIT(1)

    !$ACC END DATA
#endif

    IF (LHOOK) CALL DR_HOOK('LEINVAD_DGEMM',1,ZHOOK_HANDLE)
    !     ------------------------------------------------------------------
    END ASSOCIATE
  END SUBROUTINE LEINVAD
END MODULE LEINVAD_MOD
