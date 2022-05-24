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
SUBROUTINE LEDIR(KF_FS,KLED2,PAIA,POA1,KMODE)

!**** *LEDIR* - Direct Legendre transform.

!     Purpose.
!     --------
!        Direct Legendre tranform of state variables.

!**   Interface.
!     ----------
!        CALL LEDIR(...)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFC - number of field to transform
!                              PAIA - antisymmetric part of Fourier
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

USE PARKIND_ECTRANS  ,ONLY : JPIM ,JPIB    ,JPRB,  JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX,R_NTMAX
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_FIELDS      ,ONLY : F, &
     & ZAA,DZBST,DZCAT,ZAS,DZCST,&
     & ZAA0,DZBST0,DZCAT0,ZAS0,DZCST0,&
     & TDZAA,TDZAS,KMLOC0
USE TPM_DISTR
USE TPM_GEN, ONLY: NOUT
USE TPM_FLT
USE BUTTERFLY_ALG_MOD
USE CUDA_GEMM_BATCHED_MOD!!, ONLY: CUDA_TCGEMM_BATCHED, CUDA_GEMM_BATCHED
USE CUBLAS_MOD, ONLY : CUDA_DGEMM_BATCHED
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM) :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN)  :: KLED2
INTEGER(KIND=JPIM), INTENT(IN)  :: KMODE

REAL(KIND=JPRBT),    INTENT(IN)  :: PAIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: POA1(:,:,:)

!     LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IA, IS, ISL, J, JK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER :: ISTAT

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)

!$ACC DATA &
!$ACC& PRESENT(F,F%RW) &
!$ACC& PRESENT(D,D_NUMP,D_MYMS,R,R_NDGNH,G,G_NDGLU,R_NSMAX,R_NTMAX) &
!$ACC& PRESENT(PAIA) &
!$ACC& PRESENT(ZAA,ZAS,DZBST,DZCST,DZCAT) &
!$ACC& PRESENT(POA1,dzbst0,dzcat0,dzbst0,dzcst0) !&


! anti-symmetric

IF ( KMODE == -1 ) THEN

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISL) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,2*KF_FS
    KM = D_MYMS(KMLOC)
    IF (KM /= 0) THEN
      ISL = R_NDGNH-G_NDGLU(KM)+1

      !$ACC LOOP SEQ
      DO J=1,G_NDGLU(KM)
        DZBST((JK-1)+1+(J-1+(KMLOC-1)*R_NDGNH)*2*KF_FS)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
      ENDDO
    ENDIF
  ENDDO
END DO


! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
!$ACC HOST_DATA USE_DEVICE(ZAA,DZBST,DZCAT)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'N', &
  & 2*KF_FS, TDZAA, R_NDGNH, &
  & 1.0_JPRBT, &
  & DZBST, 2*KF_FS, R_NDGNH, &
  & ZAA, R_NDGNH, TDZAA, &
  & 0.0_JPRBT, &
  & DZCAT, 2*KF_FS, TDZAA, &
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
        POA1(JK,IA+(J-1)*2,KMLOC) = DZCAT((JK-1)+1+(J-1+(KMLOC-1)*TDZAA)*2*KF_FS)
      ENDDO
    ENDIF
  ENDDO
ENDDO

! compute m=0 in double precision:
IF(KMLOC0 > 0) THEN
   print*,'computing m=0 in double precision'

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(ISL) DEFAULT(NONE)
  DO J=1,G_NDGLU(0)
    DO JK=1,2*KF_FS,2
      ISL = R_NDGNH-G_NDGLU(0)+1
      DZBST0((JK-1)/2+1+(J-1)*2*KF_FS)=PAIA(JK,ISL+J-1,KMLOC0)*F%RW(ISL+J-1)
    ENDDO
  ENDDO


  ! Get C in transpose format to get better memory access patterns later
  !C=A*B =>
  ! C^T=B^T*A^T

  !$ACC HOST_DATA USE_DEVICE(ZAA0,DZBST0,DZCAT0)
  CALL CUDA_DGEMM_BATCHED( &
    & 'N', 'N', &
    & 2*KF_FS, TDZAA, R_NDGNH, &
    & 1.0_JPRD, &
    & DZBST0, 2*KF_FS, R_NDGNH, &
    & ZAA0, R_NDGNH, TDZAA, &
    & 0.0_JPRD, &
    & DZCAT0, 2*KF_FS, TDZAA, &
    & 1)
  !$ACC END HOST_DATA

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IA) DEFAULT(NONE)
  DO J=1,(R_NSMAX+2)/2
    DO JK=1,2*KF_FS,2
      IA  = 1+MOD(R_NTMAX+2,2)
      POA1(JK,IA+(J-1)*2,KMLOC0) = DZCAT0((JK-1)/2+1+(J-1)*2*KF_FS)
    ENDDO
  ENDDO
ENDIF

ELSE

! symmetric

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISL) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,2*KF_FS
    KM = D_MYMS(KMLOC)
    IF (KM /= 0) THEN
      ISL = R_NDGNH-G_NDGLU(KM)+1

      !$ACC LOOP SEQ
      DO J=1,G_NDGLU(KM)
        DZBST(JK+(J-1+(KMLOC-1)*R_NDGNH)*2*KF_FS)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
      ENDDO
    ENDIF
  ENDDO
END DO

! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
!$ACC HOST_DATA USE_DEVICE(ZAS,DZBST,DZCST)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'N', &
  & 2*KF_FS, TDZAS, R_NDGNH, &
  & 1.0_JPRBT, &
  & DZBST, 2*KF_FS, R_NDGNH, &
  & ZAS, R_NDGNH, TDZAS, &
  & 0.0_JPRBT, &
  & DZCST, 2*KF_FS, TDZAS, &
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
        POA1(JK,IS+(J-1)*2,KMLOC) = DZCST(JK+(J-1+(KMLOC-1)*TDZAS)*2*KF_FS)
      ENDDO
    ENDIF
  ENDDO
ENDDO

IF(KMLOC0 > 0) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(ISL) DEFAULT(NONE)
  DO J=1,G_NDGLU(0)
    DO JK=1,2*KF_FS,2
      ISL = R_NDGNH-G_NDGLU(0)+1
      DZBST0((JK-1)/2+1+(J-1)*2*KF_FS)=PAIA(JK,ISL+J-1,KMLOC0)*F%RW(ISL+J-1)
    ENDDO
  ENDDO

  ! Get C in transpose format to get better memory access patterns later
  !C=A*B =>
  ! C^T=B^T*A^T

  !$ACC host_data use_device(ZAS0,DZBST0,DZCST0)
  call CUDA_DGEMM_BATCHED( &
    & 'N', 'N', &
    & 2*KF_FS, TDZAS, R_NDGNH, &
    & 1.0_JPRD, &
    & DZBST0, 2*KF_FS, R_NDGNH, &
    & ZAS0, R_NDGNH, TDZAS, &
    & 0.0_JPRD, &
    & DZCST0, 2*KF_FS, TDZAS, &
    & 1)
  !$ACC end host_data

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IS) DEFAULT(NONE)
  DO J=1,(R_NSMAX+3)/2
    DO JK=1,2*KF_FS,2
      IS  = 1+MOD(R_NTMAX+1,2)
      POA1(JK,IS+(J-1)*2,KMLOC0) = DZCST0((JK-1)/2+1+(J-1)*2*KF_FS)
    ENDDO
  ENDDO

ENDIF

ENDIF

!$ACC END DATA


IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
