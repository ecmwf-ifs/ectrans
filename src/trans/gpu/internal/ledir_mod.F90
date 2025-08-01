#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
#if defined CUDAGPU
#define ACC_GET_HIP_STREAM ACC_GET_CUDA_STREAM
#define OPENACC_LIB OPENACC
#endif

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

MODULE LEDIR_MOD
  USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPRD, JPIB
  USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LEDIR_STRIDES, LEDIR

  INTEGER(KIND=JPIM) :: A = 8 !Alignment
CONTAINS
  SUBROUTINE LEDIR_STRIDES(KF_FS,IOUT_STRIDES0,IOUT_SIZE,IIN_STRIDES0,IIN_SIZE,&
                           IOUT0_STRIDES0,IOUT0_SIZE,IIN0_STRIDES0,IIN0_SIZE)
    USE TPM_DIM,         ONLY: R
    USE TPM_DISTR,       ONLY: D

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS

    INTEGER(KIND=JPIM), OPTIONAL :: IOUT_STRIDES0
    INTEGER(KIND=JPIB), OPTIONAL :: IOUT_SIZE
    INTEGER(KIND=JPIM), OPTIONAL :: IIN_STRIDES0
    INTEGER(KIND=JPIB), OPTIONAL :: IIN_SIZE
    INTEGER(KIND=JPIM), OPTIONAL :: IOUT0_STRIDES0, IOUT0_SIZE
    INTEGER(KIND=JPIM), OPTIONAL :: IIN0_STRIDES0, IIN0_SIZE

    ASSOCIATE(D_OFFSETS_GEMM1=>D%OFFSETS_GEMM1, D_OFFSETS_GEMM2=>D%OFFSETS_GEMM2)

    IF (PRESENT(IOUT_STRIDES0)) &
        IOUT_STRIDES0 = ALIGN(2*KF_FS,A)
    IF (PRESENT(IOUT_SIZE)) &
        IOUT_SIZE = IOUT_STRIDES0*D_OFFSETS_GEMM2(D%NUMP+1)
    IF (PRESENT(IIN_STRIDES0)) &
        IIN_STRIDES0 = ALIGN(2*KF_FS,A)
    IF (PRESENT(IIN_SIZE)) &
        IIN_SIZE = IIN_STRIDES0*D_OFFSETS_GEMM1(D%NUMP+1)
    IF (PRESENT(IOUT0_STRIDES0)) &
        IOUT0_STRIDES0 = ALIGN(KF_FS,A)
    IF (PRESENT(IOUT0_SIZE)) &
        IOUT0_SIZE = IOUT0_STRIDES0 * ALIGN(MAX((R%NTMAX+2)/2,(R%NTMAX+3)/2),A)
    IF (PRESENT(IIN0_STRIDES0)) &
        IIN0_STRIDES0 = ALIGN(KF_FS,A)
    IF (PRESENT(IIN0_SIZE)) &
        IIN0_SIZE = IIN0_STRIDES0 * ALIGN(R%NDGNH,A)

    END ASSOCIATE
  END SUBROUTINE

  SUBROUTINE LEDIR(ALLOCATOR,ZINPS,ZINPA,ZINPS0,ZINPA0,ZOUT,ZOUT0,POA1,KF_FS)
    !**** *LEDIR* - Direct Legendre transform.

    !     Purpose.
    !     --------
    !        Direct Legendre tranform of state variables.

    !**   Interface.
    !     ----------
    !        CALL LEDIR(...)

    !        Explicit arguments :  KM - zonal wavenumber
    !        --------------------  KFC - number of field to transform
    !                              fields for zonal wavenumber KM
    !                              PSIA - symmetric part of Fourier
    !                              fields for zonal wavenumber KM
    !                              POA1 -  spectral
    !                              fields for zonal wavenumber KM

    !        Implicit arguments :  None.
    !        --------------------

    !     Method.
    !     -------   use butterfly or dgemm

    !     Externals.   
    !     ----------

    !     Reference.
    !     ----------
    !        ECMWF Research Department documentation of the IFS

    !     Author.
    !     -------
    !          Nils Wedi + Mats Hamrud + George Modzynski

    !     Modifications.
    !     --------------
    !        J.Hague : Oct 2012 DR_HOOK round calls to DGEMM:
    !      F. Vana  05-Mar-2015  Support for single precision
    !     ------------------------------------------------------------------

    USE TPM_GEN,                     ONLY: LSYNC_TRANS, NOUT, NCUR_RESOL
    USE YOMHOOK,                     ONLY: LHOOK,   DR_HOOK, JPHOOK
    USE TPM_DIM,                     ONLY: R
    USE TPM_GEOMETRY,                ONLY: G
    USE TPM_FIELDS_GPU,              ONLY: FG
    USE TPM_DISTR,                   ONLY: D
    USE HICBLAS_MOD,                 ONLY: HIP_DGEMM_BATCHED, &
      &                                    HIP_DGEMM_GROUPED, HIP_SGEMM_GROUPED
    USE MPL_MODULE,                  ONLY: MPL_BARRIER,MPL_ALL_MS_COMM
    USE TPM_STATS,                   ONLY: GSTATS => GSTATS_NVTX
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_LONG, C_LOC
#ifdef ACCGPU
    USE OPENACC_LIB, ONLY: ACC_GET_HIP_STREAM
#endif

#ifdef TRANS_SINGLE
#define HIP_GEMM HIP_SGEMM_GROUPED
#else
#define HIP_GEMM HIP_DGEMM_GROUPED
#endif

    IMPLICIT NONE

    !     DUMMY ARGUMENTS
    REAL(KIND=JPRBT), INTENT(IN) :: ZINPS(:), ZINPA(:)
    REAL(KIND=JPRD), INTENT(IN) :: ZINPS0(:), ZINPA0(:)
    REAL(KIND=JPRBT), INTENT(INOUT) :: ZOUT(:)
    REAL(KIND=JPRD), INTENT(INOUT) ::  ZOUT0(:)
    REAL(KIND=JPRBT),  INTENT(OUT), POINTER :: POA1(:,:,:)
    INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR

    !     LOCAL VARIABLES
    INTEGER(KIND=JPIM)  :: KM
    INTEGER(KIND=JPIM)  :: KMLOC
    INTEGER(KIND=JPIM) :: IA, IS, J, IMLOC0(1)
    INTEGER(KIND=JPIM)  :: KS(D%NUMP), NS(D%NUMP)
    INTEGER(KIND=JPIB)  :: AOFFSETS(D%NUMP), BOFFSETS(D%NUMP), COFFSETS(D%NUMP)
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    INTEGER(KIND=JPIM) :: JF

    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0
    INTEGER(KIND=JPIB)  :: IOUT_STRIDES1
    INTEGER(KIND=JPIM)  :: IIN_STRIDES0
    INTEGER(KIND=JPIB)  :: IIN_STRIDES1
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_STRIDES1
    INTEGER(KIND=JPIM)  :: IIN0_STRIDES0, IIN0_STRIDES1

    INTEGER(KIND=C_LONG) :: HIP_STREAM

    ASSOCIATE(D_NUMP=>D%NUMP, R_NSMAX=>R%NSMAX, R_NTMAX=>R%NTMAX, G_NDGLU=>G%NDGLU, &
            & D_MYMS=>D%MYMS, D_OFFSETS_GEMM1=>D%OFFSETS_GEMM1, &
            & D_OFFSETS_GEMM2=>D%OFFSETS_GEMM2, &
            & ZAA=>FG%ZAA, ZAS=>FG%ZAS, ZAA0=>FG%ZAA0, ZAS0=>FG%ZAS0)
    IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)

#ifdef ACCGPU
    HIP_STREAM = INT(ACC_GET_HIP_STREAM(1_C_INT), C_LONG)
#endif
#ifdef OMPGPU
    HIP_STREAM = 0_C_LONG
#endif

    CALL LEDIR_STRIDES(KF_FS,IOUT_STRIDES0,IOUT_STRIDES1,IIN_STRIDES0,IIN_STRIDES1,&
                       IOUT0_STRIDES0,IOUT0_STRIDES1,IIN0_STRIDES0,IIN0_STRIDES1)

#ifdef OMPGPU
    !$OMP TARGET DATA &
    !$OMP& MAP(PRESENT,ALLOC:ZINPS,ZINPA,ZOUT,ZINPS0,ZINPA0,ZOUT0) &
    !$OMP& MAP(PRESENT,ALLOC:D,D_MYMS,D_NUMP,R,R_NTMAX,R_NSMAX) &
    !$OMP& MAP(PRESENT,ALLOC:ZAA,ZAS,POA1,D_OFFSETS_GEMM1,D_OFFSETS_GEMM2)
#endif
#ifdef ACCGPU
    !$ACC DATA &
    !$ACC& PRESENT(ZINPS,ZINPA,ZOUT,ZINPS0,ZINPA0,ZOUT0) &
    !$ACC& PRESENT(D,D_MYMS,D_NUMP,R,R_NTMAX,R_NSMAX) &
    !$ACC& PRESENT(ZAA,ZAS,POA1,D_OFFSETS_GEMM1,D_OFFSETS_GEMM2)
#endif

    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(430,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(430,1)
    ENDIF
    CALL GSTATS(414,0)

    ! anti-symmetric
    IMLOC0 = FINDLOC(D_MYMS,0)
    IF(IMLOC0(1) > 0) THEN
      ! compute m=0 in double precision:
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAA0,ZINPA0,ZOUT0)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAA0,ZINPA0,ZOUT0)
#endif
      CALL HIP_DGEMM_BATCHED( &
        & 'N', 'N', &
        & KF_FS, (R_NSMAX+2)/2, G_NDGLU(0), &
        & 1.0_JPRD, &
        & ZINPA0, IIN0_STRIDES0, 0, &
        & ZAA0, SIZE(ZAA0,1), 0, &
        & 0.0_JPRD, &
        & ZOUT0, IOUT0_STRIDES0, 0, &
        & 1, HIP_STREAM, C_LOC(ALLOCATOR%PTR))
#ifdef OMPGPU
    !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
    !$ACC END HOST_DATA
#endif
    ENDIF
    ! Get C in transpose format to get better memory access patterns later
    !C=A*B =>
    ! C^T=B^T*A^T
    DO KMLOC=1,D_NUMP
      KM = D_MYMS(KMLOC)
      NS(KMLOC) = (R_NSMAX-KM+2)/2
      KS(KMLOC) = G_NDGLU(KM)
      AOFFSETS(KMLOC) = IIN_STRIDES0*D_OFFSETS_GEMM1(KMLOC)
      BOFFSETS(KMLOC) = D%OFFSETS_GEMM_MATRIX(KMLOC)
      COFFSETS(KMLOC) = IOUT_STRIDES0*D_OFFSETS_GEMM2(KMLOC)
    ENDDO
    IF(IMLOC0(1) > 0) THEN
      NS(IMLOC0(1)) = 0
      KS(IMLOC0(1)) = 0
    ENDIF
#ifdef OMPGPU
    !$OMP TARGET DATA USE_DEVICE_PTR(ZAA,ZINPA,ZOUT)
#endif
#ifdef ACCGPU
    !$ACC HOST_DATA USE_DEVICE(ZAA,ZINPA,ZOUT)
#endif
    CALL HIP_GEMM( &
      & NCUR_RESOL, 21, & ! unique identifier
      & 'N', 'N', &
      & 2*KF_FS, NS(:), KS(:), &
      & 1.0_JPRBT, &
      & ZINPA, IIN_STRIDES0, AOFFSETS, &
      & ZAA, D%LEGENDRE_MATRIX_STRIDES, BOFFSETS, &
      & 0.0_JPRBT, &
      & ZOUT, IOUT_STRIDES0, COFFSETS, &
      & D_NUMP, HIP_STREAM, C_LOC(ALLOCATOR%PTR))
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
      CALL GSTATS(434,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(434,1)
    ENDIF
    CALL GSTATS(414,1)

#ifdef OMPGPU
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) PRIVATE(KM,IA) &
    !$OMP& SHARED(D,R,KF_FS,IOUT_STRIDES0,ZOUT,IOUT0_STRIDES0,ZOUT0,POA1) &
    !$OMP& MAP(TO:KF_FS,IOUT_STRIDES0)
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IA,J) FIRSTPRIVATE(KF_FS,IOUT_STRIDES0,IOUT0_STRIDES0) DEFAULT(NONE) &
#ifndef _CRAYFTN
    !$ACC& ASYNC(1)
#else
    !$ACC&
#endif
#endif
    DO KMLOC=1,D_NUMP
      DO JF=1,2*KF_FS
        KM = D_MYMS(KMLOC)
        IA  = 1+MOD(R_NTMAX-KM+2,2)
        IF (KM /= 0) THEN
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX-KM+2)/2
            POA1(JF,IA+1+(J-1)*2,KMLOC) = ZOUT(JF+(J-1)*IOUT_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IOUT_STRIDES0)
          ENDDO
        ELSEIF (MOD(JF-1,2) == 0) THEN
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX+2)/2
            POA1(JF,IA+1+(J-1)*2,KMLOC) = ZOUT0((JF-1)/2+1+(J-1)*IOUT0_STRIDES0)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    ! symmetric

    IF (LSYNC_TRANS) THEN
#ifdef ACCGPU
      !$ACC WAIT(1)
#endif
      CALL GSTATS(430,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(430,1)
    ENDIF
    CALL GSTATS(414,0)

    IF(IMLOC0(1) > 0) THEN
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAS0,ZINPS0,ZOUT0)
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAS0,ZINPS0,ZOUT0)
#endif
      ! compute m=0 in double precision:
      call HIP_DGEMM_BATCHED( &
        & 'N', 'N', &
        & KF_FS, (R_NSMAX+3)/2, G_NDGLU(0), &
        & 1.0_JPRD, &
        & ZINPS0, IIN0_STRIDES0, 0, &
        & ZAS0, SIZE(ZAS0,1), 0, &
        & 0.0_JPRD, &
        & ZOUT0, IOUT0_STRIDES0, 0, &
        & 1, HIP_STREAM, C_LOC(ALLOCATOR%PTR))
#ifdef OMPGPU
    !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
    !$ACC END HOST_DATA
#endif
    ENDIF

    ! Get C in transpose format to get better memory access patterns later
    !C=A*B =>
    ! C^T=B^T*A^T
    DO KMLOC=1,D_NUMP
      KM = D_MYMS(KMLOC)
      NS(KMLOC) = (R_NSMAX-KM+3)/2
      KS(KMLOC) = G_NDGLU(KM)
      AOFFSETS(KMLOC) = IIN_STRIDES0*D_OFFSETS_GEMM1(KMLOC)
      BOFFSETS(KMLOC) = D%OFFSETS_GEMM_MATRIX(KMLOC)
      COFFSETS(KMLOC) = IOUT_STRIDES0*D_OFFSETS_GEMM2(KMLOC)
    ENDDO
    IF(IMLOC0(1) > 0) THEN
      NS(IMLOC0(1)) = 0
      KS(IMLOC0(1)) = 0
    ENDIF
#ifdef OMPGPU
    !$OMP TARGET DATA USE_DEVICE_PTR(ZAS,ZINPS,ZOUT)
#endif
#ifdef ACCGPU
    !$ACC HOST_DATA USE_DEVICE(ZAS,ZINPS,ZOUT)
#endif
    CALL HIP_GEMM( &
      & NCUR_RESOL, 22, & ! unique identifier
      & 'N', 'N', &
      & 2*KF_FS, NS(:), KS(:), &
      & 1.0_JPRBT, &
      & ZINPS, IIN_STRIDES0, AOFFSETS, &
      & ZAS, D%LEGENDRE_MATRIX_STRIDES, BOFFSETS, &
      & 0.0_JPRBT, &
      & ZOUT, IOUT_STRIDES0, COFFSETS, &
      & D_NUMP, HIP_STREAM, C_LOC(ALLOCATOR%PTR))
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
      CALL GSTATS(434,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(434,1)
    ENDIF
    CALL GSTATS(414,1)

#ifdef OMPGPU
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(KM,IS) &
    !$OMP& SHARED(D,R,KF_FS,IOUT_STRIDES0,ZOUT,POA1)
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IS) FIRSTPRIVATE(KF_FS,IOUT_STRIDES0,IOUT0_STRIDES0) &
    !$ACC& DEFAULT(NONE) &
#ifndef _CRAYFTN
    !$ACC& ASYNC(1)
#else
    !$ACC&
#endif
#endif
    DO KMLOC=1,D_NUMP
      DO JF=1,2*KF_FS
        KM = D_MYMS(KMLOC)
        IS  = 1+MOD(R_NTMAX-KM+1,2)
        IF (KM /= 0) THEN
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX-KM+3)/2
            POA1(JF,IS+1+(J-1)*2,KMLOC) = ZOUT(JF+(J-1)*IOUT_STRIDES0+D_OFFSETS_GEMM2(KMLOC)*IOUT_STRIDES0)
          ENDDO
        ELSEIF (MOD(JF-1,2) == 0) THEN
#ifdef ACCGPU
          !$ACC LOOP SEQ
#endif
          DO J=1,(R_NSMAX+3)/2
            POA1(JF,IS+1+(J-1)*2,KMLOC) = ZOUT0((JF-1)/2+1+(J-1)*IOUT0_STRIDES0)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
#ifdef OMPGPU
    !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
    !$ACC WAIT(1)

    !$ACC END DATA
#endif

    IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
    !     ------------------------------------------------------------------
    END ASSOCIATE
  END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
