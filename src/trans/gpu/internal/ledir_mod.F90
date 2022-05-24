! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEDIR_MOD
CONTAINS
SUBROUTINE LEDIR(FOUBUF,POA1,KF_FS,KF_UV)

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

USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE PARKIND_ECTRANS  ,ONLY : JPIM, JPRBT, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX,R_NTMAX,R_NDGL
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_FIELDS      ,ONLY : F,ZAA,ZAS,ZAA0,ZAS0,TDZAA,TDZAS,KMLOC0
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS, D_NPNTGTB1
USE CUDA_GEMM_BATCHED_MOD
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC


IMPLICIT NONE


!     DUMMY ARGUMENTS
REAL(KIND=JPRBT), ALLOCATABLE, INTENT(INOUT) :: FOUBUF(:)
REAL(KIND=JPRBT),  INTENT(OUT) :: POA1(:,:,:)
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS, KF_UV

!     LOCAL VARIABLES
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM) :: IA, IS, ISL, J, JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRBT) :: PAIA, PAIS

INTEGER(KIND=JPIM) :: IGLS,  JF, JGL
INTEGER(KIND=JPIM) :: OFFSET1, OFFSET2
REAL(KIND=JPRBT), ALLOCATABLE :: ZINPS(:), ZINPA(:), ZOUT(:)
REAL(KIND=JPRD), ALLOCATABLE :: ZINP0(:), ZOUT0(:)

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)

ALLOCATE(ZINPA(2*KF_FS*R_NDGNH*D_NUMP))
ALLOCATE(ZINPS(2*KF_FS*R_NDGNH*D_NUMP))
ALLOCATE(ZOUT(2*KF_FS*TDZAS*D_NUMP))
ALLOCATE(ZINP0(KF_FS*R_NDGNH))
ALLOCATE(ZOUT0(KF_FS*TDZAS))

!$ACC DATA &
!$ACC& CREATE(ZINPS,ZINPA,ZOUT,ZINP0,ZOUT0) &
!$ACC& PRESENT(F,F%RW) &
!$ACC& PRESENT(D,D_MYMS,R,G,G_NDGLU) &
!$ACC& PRESENT(ZAA,ZAS,POA1) &
!$ACC& PRESENT(D_NPNTGTB1)

! TODO this doesn't make sense that we need it (???)
!$ACC KERNELS
ZINPS(:) = 0
ZINPA(:) = 0
!$ACC END KERNELS

!$ACC DATA PRESENT(FOUBUF)
!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(KM,ISL,IGLS,OFFSET1,OFFSET2,JGL,PAIA,PAIS)
DO KMLOC=1,D_NUMP
  DO JGL=1,R_NDGNH
    DO JF=1,KF_FS*2
      KM = D_MYMS(KMLOC)
      ISL = R_NDGNH-G_NDGLU(KM)+1
      IF (JGL >= ISL) THEN
        !(DO JGL=ISL,R_NDGNH)
        IGLS = R_NDGL+1-JGL
        OFFSET1 = D_NPNTGTB1(KMLOC,JGL )*2*KF_FS
        OFFSET2 = D_NPNTGTB1(KMLOC,IGLS)*2*KF_FS
        PAIA = FOUBUF(OFFSET1+JF)-FOUBUF(OFFSET2+JF)
        PAIS = FOUBUF(OFFSET1+JF)+FOUBUF(OFFSET2+JF)
        IF (JF <= 4*KF_UV) THEN
            ! Multiply in case of velocity
          PAIA = PAIA*F%RACTHE(JGL)
          PAIS = PAIS*F%RACTHE(JGL)
        ENDIF
        ZINPA(JF+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_FS)=PAIA*F%RW(JGL)
        ZINPS(JF+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_FS)=PAIS*F%RW(JGL)
      ENDIF
    ENDDO
  ENDDO
END DO
!$ACC END DATA

!$ACC EXIT DATA DELETE(FOUBUF)
DEALLOCATE(FOUBUF)

! anti-symmetric

IF (LSYNC_TRANS) THEN
  CALL GSTATS(430,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(430,1)
ENDIF
CALL GSTATS(414,0)
! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
DO KMLOC=1,D_NUMP
  CALL CUDA_GEMM_BATCHED( &
    & CUBLAS_OP_N, CUBLAS_OP_N, &
    & 2*KF_FS, TDZAA, R_NDGNH, &
    & 1.0_JPRBT, &
    & ZINPA((KMLOC-1)*R_NDGNH*2*KF_FS+1:), 2*KF_FS, R_NDGNH*2*KF_FS, &
    & ZAA(:,:,KMLOC), R_NDGNH, TDZAA*R_NDGNH, &
    & 0.0_JPRBT, &
    & ZOUT((KMLOC-1)*TDZAA*2*KF_FS+1:), 2*KF_FS, TDZAA*2*KF_FS, &
    & 1)
ENDDO
IF (LSYNC_TRANS) THEN
  CALL GSTATS(434,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(434,1)
ENDIF
CALL GSTATS(414,1)

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IA) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,2*KF_FS
    KM = D_MYMS(KMLOC)
    IF (KM /= 0) THEN
      IA  = 1+MOD(R_NTMAX-KM+2,2)
      !$ACC LOOP SEQ
      DO J=1,(R%NSMAX-KM+2)/2
        POA1(JK,IA+1+(J-1)*2,KMLOC) = ZOUT(JK+(J-1+(KMLOC-1)*TDZAA)*2*KF_FS)
      ENDDO
    ENDIF
  ENDDO
ENDDO

! compute m=0 in double precision:
IF(KMLOC0 > 0) THEN
   print*,'computing m=0 in double precision'

  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE)
  DO JGL=1,G_NDGLU(0)
    DO JF=1,2*KF_FS,2
      ZINP0((JF-1)/2+1+(JGL-1)*KF_FS) &
          & = ZINPA((JF-1)+1+(JGL-1+(KMLOC0-1)*R_NDGNH)*2*KF_FS)
    ENDDO
  ENDDO


  ! Get C in transpose format to get better memory access patterns later
  !C=A*B =>
  ! C^T=B^T*A^T

  CALL CUDA_GEMM_BATCHED( &
    & CUBLAS_OP_N, CUBLAS_OP_N, &
    & KF_FS, TDZAA, R_NDGNH, &
    & 1.0_JPRD, &
    & ZINP0, KF_FS, R_NDGNH*KF_FS, &
    & ZAA0, R_NDGNH, TDZAA*R_NDGNH, &
    & 0.0_JPRD, &
    & ZOUT0, KF_FS, TDZAA*KF_FS, &
    & 1)

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IA) DEFAULT(NONE)
  DO J=1,(R_NSMAX+2)/2
    DO JK=1,2*KF_FS,2
      IA  = 1+MOD(R_NTMAX+2,2)
      POA1(JK,IA+1+(J-1)*2,KMLOC0) = ZOUT0((JK-1)/2+1+(J-1)*KF_FS)
    ENDDO
  ENDDO
ENDIF

! symmetric

IF (LSYNC_TRANS) THEN
  CALL GSTATS(430,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(430,1)
ENDIF
CALL GSTATS(414,0)
! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
DO KMLOC=1,D_NUMP
  CALL CUDA_GEMM_BATCHED( &
    & CUBLAS_OP_N, CUBLAS_OP_N, &
    & 2*KF_FS, TDZAS, R_NDGNH, &
    & 1.0_JPRBT, &
    & ZINPS((KMLOC-1)*R_NDGNH*2*KF_FS+1:), 2*KF_FS, R_NDGNH*2*KF_FS, &
    & ZAS(:,:,KMLOC), R_NDGNH, TDZAS*R_NDGNH, &
    & 0.0_JPRBT, &
    & ZOUT((KMLOC-1)*TDZAS*2*KF_FS+1:), 2*KF_FS, TDZAS*2*KF_FS, &
    & 1)
ENDDO
IF (LSYNC_TRANS) THEN
  CALL GSTATS(434,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(434,1)
ENDIF
CALL GSTATS(414,1)

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IS) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,2*KF_FS
    KM = D_MYMS(KMLOC)
    IF (KM /= 0) THEN
      IS  = 1+MOD(R_NTMAX-KM+1,2)
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX-KM+3)/2
        POA1(JK,IS+1+(J-1)*2,KMLOC) = ZOUT(JK+(J-1+(KMLOC-1)*TDZAS)*2*KF_FS)
      ENDDO
    ENDIF
  ENDDO
ENDDO

IF(KMLOC0 > 0) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE)
  DO JGL=1,G_NDGLU(0)
    DO JF=1,2*KF_FS,2
      ZINP0((JF-1)/2+1+(JGL-1)*KF_FS) &
          & = ZINPS((JF-1)+1+(JGL-1+(KMLOC0-1)*R_NDGNH)*2*KF_FS)
    ENDDO
  ENDDO

  ! Get C in transpose format to get better memory access patterns later
  !C=A*B =>
  ! C^T=B^T*A^T

  call CUDA_GEMM_BATCHED( &
    & CUBLAS_OP_N, CUBLAS_OP_N, &
    & KF_FS, TDZAS, R_NDGNH, &
    & 1.0_JPRD, &
    & ZINP0, KF_FS, R_NDGNH*KF_FS, &
    & ZAS0, R_NDGNH, TDZAS*R_NDGNH, &
    & 0.0_JPRD, &
    & ZOUT0, KF_FS, TDZAS*KF_FS, &
    & 1)

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IS) DEFAULT(NONE)
  DO J=1,(R_NSMAX+3)/2
    DO JK=1,2*KF_FS,2
      IS  = 1+MOD(R_NTMAX+1,2)
      POA1(JK,IS+1+(J-1)*2,KMLOC0) = ZOUT0((JK-1)/2+1+(J-1)*KF_FS)
    ENDDO
  ENDDO

ENDIF

!$ACC END DATA
DEALLOCATE(ZINPA)
DEALLOCATE(ZINPS)
DEALLOCATE(ZOUT)
DEALLOCATE(ZINP0)
DEALLOCATE(ZOUT0)


IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
