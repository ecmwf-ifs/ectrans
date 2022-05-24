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
SUBROUTINE LEDIR(KF_FS,KF_UV,POA1)

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
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS, D_NPROCL,D_NPNTGTB1,D_NSTAGT1B
USE TPM_TRANS       ,ONLY : FOUBUF
USE CUDA_GEMM_BATCHED_MOD!!, ONLY: CUDA_TCGEMM_BATCHED, CUDA_GEMM_BATCHED
USE CUBLAS_MOD, ONLY : CUDA_DGEMM_BATCHED
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC


IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS, KF_UV
REAL(KIND=JPRBT),    INTENT(OUT) :: POA1(:,:,:)

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

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEDIR BARRIER')
ENDIF
CALL GSTATS(452,0)

ALLOCATE(ZINPA(2*KF_FS*R_NDGNH*D_NUMP))
ALLOCATE(ZINPS(2*KF_FS*R_NDGNH*D_NUMP))
ALLOCATE(ZOUT(2*KF_FS*TDZAS*D_NUMP))
ALLOCATE(ZINP0(2*KF_FS*R_NDGNH))
ALLOCATE(ZOUT0(2*KF_FS*TDZAS))

!$ACC DATA &
!$ACC& CREATE(ZINPS,ZINPA,ZOUT,ZINP0,ZOUT0) &
!$ACC& PRESENT(F,F%RW) &
!$ACC& PRESENT(D,D_MYMS,R,G,G_NDGLU) &
!$ACC& PRESENT(ZAA,ZAS,POA1) &
!$ACC& PRESENT(FOUBUF,D_NPNTGTB1,D_NSTAGT1B,D_NPROCL)

! TODO this doesn't make sense that we need it (???)
!$ACC KERNELS
ZINPS(:) = 0
ZINPA(:) = 0
!$ACC END KERNELS

!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(2) PRIVATE(KM,ISL,IGLS,OFFSET1,OFFSET2,JGL,PAIA,PAIS)
DO KMLOC=1,D_NUMP
  DO JF=1,KF_FS*2
    KM = D_MYMS(KMLOC)
    ISL = R_NDGNH-G_NDGLU(KM)+1
    !$ACC LOOP SEQ
    DO JGL=ISL,R_NDGNH
      IGLS = R_NDGL+1-JGL
      OFFSET1 = (D_NSTAGT1B(D_NPROCL(JGL) )+D_NPNTGTB1(KMLOC,JGL ))*2*KF_FS
      OFFSET2 = (D_NSTAGT1B(D_NPROCL(IGLS))+D_NPNTGTB1(KMLOC,IGLS))*2*KF_FS
      PAIA = FOUBUF(OFFSET1+JF)-FOUBUF(OFFSET2+JF)
      PAIS = FOUBUF(OFFSET1+JF)+FOUBUF(OFFSET2+JF)
      IF (JF .LE. 4*KF_UV) THEN
          ! Multiply in case of velocity
        PAIA = PAIA*F%RACTHE(JGL)
        PAIS = PAIS*F%RACTHE(JGL)
      ENDIF
      ZINPA(JF+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_FS)=PAIA*F%RW(JGL)
      ZINPS(JF+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_FS)=PAIS*F%RW(JGL)
    ENDDO
  ENDDO
END DO

! anti-symmetric

! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
!$ACC HOST_DATA USE_DEVICE(ZAA,ZINPA,ZOUT)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'N', &
  & 2*KF_FS, TDZAA, R_NDGNH, &
  & 1.0_JPRBT, &
  & ZINPA, 2*KF_FS, R_NDGNH, &
  & ZAA, R_NDGNH, TDZAA, &
  & 0.0_JPRBT, &
  & ZOUT, 2*KF_FS, TDZAA, &
  & D_NUMP)
!$ACC END HOST_DATA

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
      ZINP0((JF-1)/2+1+(JGL-1)*2*KF_FS) &
          & = ZINPA((JF-1)+1+(JGL-1+(KMLOC0-1)*R_NDGNH)*2*KF_FS)
    ENDDO
  ENDDO


  ! Get C in transpose format to get better memory access patterns later
  !C=A*B =>
  ! C^T=B^T*A^T

  !$ACC HOST_DATA USE_DEVICE(ZAA0,ZINP0,ZOUT0)
  CALL CUDA_DGEMM_BATCHED( &
    & 'N', 'N', &
    & 2*KF_FS, TDZAA, R_NDGNH, &
    & 1.0_JPRD, &
    & ZINP0, 2*KF_FS, R_NDGNH, &
    & ZAA0, R_NDGNH, TDZAA, &
    & 0.0_JPRD, &
    & ZOUT0, 2*KF_FS, TDZAA, &
    & 1)
  !$ACC END HOST_DATA

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IA) DEFAULT(NONE)
  DO J=1,(R_NSMAX+2)/2
    DO JK=1,2*KF_FS,2
      IA  = 1+MOD(R_NTMAX+2,2)
      POA1(JK,IA+1+(J-1)*2,KMLOC0) = ZOUT0((JK-1)/2+1+(J-1)*2*KF_FS)
    ENDDO
  ENDDO
ENDIF

! symmetric

! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
!$ACC HOST_DATA USE_DEVICE(ZAS,ZINPS,ZOUT)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'N', &
  & 2*KF_FS, TDZAS, R_NDGNH, &
  & 1.0_JPRBT, &
  & ZINPS, 2*KF_FS, R_NDGNH, &
  & ZAS, R_NDGNH, TDZAS, &
  & 0.0_JPRBT, &
  & ZOUT, 2*KF_FS, TDZAS, &
  & D_NUMP)
!$ACC END HOST_DATA

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
      ZINP0((JF-1)/2+1+(JGL-1)*2*KF_FS) &
          & = ZINPS((JF-1)+1+(JGL-1+(KMLOC0-1)*R_NDGNH)*2*KF_FS)
    ENDDO
  ENDDO

  ! Get C in transpose format to get better memory access patterns later
  !C=A*B =>
  ! C^T=B^T*A^T

  !$ACC host_data use_device(ZAS0,ZINP0,ZOUT0)
  call CUDA_DGEMM_BATCHED( &
    & 'N', 'N', &
    & 2*KF_FS, TDZAS, R_NDGNH, &
    & 1.0_JPRD, &
    & ZINP0, 2*KF_FS, R_NDGNH, &
    & ZAS0, R_NDGNH, TDZAS, &
    & 0.0_JPRD, &
    & ZOUT0, 2*KF_FS, TDZAS, &
    & 1)
  !$ACC end host_data

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IS) DEFAULT(NONE)
  DO J=1,(R_NSMAX+3)/2
    DO JK=1,2*KF_FS,2
      IS  = 1+MOD(R_NTMAX+1,2)
      POA1(JK,IS+1+(J-1)*2,KMLOC0) = ZOUT0((JK-1)/2+1+(J-1)*2*KF_FS)
    ENDDO
  ENDDO

ENDIF

!$ACC END DATA
DEALLOCATE(ZINPA)
DEALLOCATE(ZINPS)
DEALLOCATE(ZOUT)
DEALLOCATE(ZINP0)
DEALLOCATE(ZOUT0)


IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEDIR BARRIER')
ENDIF
CALL GSTATS(452,1)
IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
