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
CONTAINS
SUBROUTINE LEINV(KFC,KF_OUT_LT,PIA,PAOA1,PSOA1)

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
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)

!        Implicit arguments :  None.
!        --------------------

!     Method.    use butterfly or dgemm
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
USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_FIELDS      ,ONLY : F, ZIA, &
     & ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
     & IZBS,ILDZBA,ILDZBS,ITDZBA,ITDZBS,&
     & IZCS,IZCST,ILDZCA,ILDZCS,ITDZCA,ITDZCS,&
     & TDZAS, IF_FS_INV, ZAMAX, ZSMAX 
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_GEN, ONLY: NOUT
USE TPM_FLT
USE BUTTERFLY_ALG_MOD
USE CUDA_GEMM_BATCHED_MOD
USE, INTRINSIC :: ISO_C_BINDING
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PSOA1(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PAOA1(:,:,:)

!     LOCAL
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, J1, IF, JGL,JK, J,JI, IRET
REAL(KIND=JPRBT), ALLOCATABLE :: ZZBS(:,:,:), ZZCSTS(:,:,:)
REAL(KIND=JPRBT), ALLOCATABLE :: ZZBA(:,:,:), ZZCSTA(:,:,:)

INTEGER(KIND=JPIM) :: ISTAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  
!*       1.1      PREPARATIONS.
IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.
IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEINV BARRIER')
ENDIF
CALL GSTATS(453,0)

ALLOCATE(ZZBA(KFC,TDZAA,D_NUMP))
ALLOCATE(ZZCSTA(KFC,R_NDGNH,D_NUMP))
ALLOCATE(ZZBS(KFC,TDZAS,D_NUMP))
ALLOCATE(ZZCSTS(KFC,R_NDGNH,D_NUMP))
!$ACC DATA CREATE(ZZBA,ZZCSTA,ZZBS,ZZCSTS)

!$ACC DATA COPYIN (D,D_MYMS,R,G,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX) &
!$ACC&     PRESENT (ZAA,ZAS) &
!$ACC&     PRESENT (ZZBA,ZZBS,ZZCSTA,ZZCSTS) &
!$ACC&     PRESENT (PIA) &
!$ACC&     PRESENT (PSOA1,PAOA1,IZBS)

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
   DO JGL=1,R_NDGNH
      DO J1=2,KFC,2
         KM = D_MYMS(KMLOC)
         IF(KM == 0)THEN
            PSOA1(J1,JGL,KMLOC) = 0.0_JPRBT
            PAOA1(J1,JGL,KMLOC) = 0.0_JPRBT
         END IF
      ENDDO
   ENDDO
   !end loop over wavenumber
END DO

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IA) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,KFC
    KM =  D_MYMS(KMLOC)
    IF(KM /= 0)THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX-KM+2)/2
        IA  = 1+MOD(R_NSMAX-KM+2,2)
        ZZBA(JK,J,KMLOC)=PIA(JK,IA+1+(J-1)*2,KMLOC)
      ENDDO
    ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX+2)/2
        IA  = 1+MOD(R_NSMAX+2,2)
        ZZBA((JK-1)/2+1,J,KMLOC)=PIA(JK,IA+1+(J-1)*2,KMLOC)
      ENDDO
    ENDIF
  ENDDO
ENDDO

! operate on full arrays, where non-relavent entries have been set to zero
! call CUDA_DGEMM_BATCHED('N','N',LDZAA,TDZBA,TDZAA,1.0_JPRB,ZAA,LDZAA,TDZAA,ZBA,LDZBA,TDZBA,0._JPRB,ZCA,LDZCA,TDZCA,D_NUMP)
! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T


! OVERLOADED FOR SINGLE AND DOUBLE PRECISION
!$ACC HOST_DATA USE_DEVICE(ZAA,ZZBA,ZZCSTA)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'T', &
  & KFC, R_NDGNH, TDZAA, &
  & 1.0_JPRBT, &
  & ZZBA, KFC, TDZAA,&
  & ZAA, R_NDGNH, TDZAA, &
  & 0._JPRBT, &
  & ZZCSTA, KFC, R_NDGNH, &
  & D_NUMP)
!$ACC END HOST_DATA

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISL) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,KFC
    KM = D_MYMS(KMLOC)
    ISL = R_NDGNH-G_NDGLU(KM)+1
    !$ACC LOOP SEQ
    DO JGL=ISL,R_NDGNH
      IF(KM /= 0) THEN
        PAOA1(JK,JGL,KMLOC) = ZZCSTA((JK-1)+1,JGL-ISL+1,KMLOC)
      ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
        PAOA1(JK,JGL,KMLOC) = ZZCSTA((JK-1)/2+1,JGL-ISL+1,KMLOC)
      ENDIF
    ENDDO
  ENDDO
ENDDO

! 2. +++++++++++++ symmetric

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IS) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,KFC
    KM =  D_MYMS(KMLOC)
    IF(KM /= 0)THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX-KM+3)/2
        IS  = 1+MOD(R_NSMAX-KM+1,2)
        ZZBS(JK,J,KMLOC)=PIA(JK,IS+1+(J-1)*2,KMLOC)
      ENDDO
    ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX+3)/2
        IS  = 1+MOD(R_NSMAX+1,2)
        ZZBS((JK-1)/2+1,J,KMLOC)=PIA(JK,IS+1+(J-1)*2,KMLOC)
      ENDDO
    ENDIF
  ENDDO
ENDDO

!C=A*B =>
! C^T=B^T*A^T

!$ACC HOST_DATA USE_DEVICE(ZAS,ZZBS,ZZCSTS)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'T', &
  & KFC, R_NDGNH, TDZAS, &
  & 1.0_JPRBT, &
  & ZZBS, KFC, TDZAS, &
  & ZAS, R_NDGNH, TDZAS, &
  & 0._JPRBT, &
  & ZZCSTS, KFC, R_NDGNH, &
  & D_NUMP)
!$ACC END HOST_DATA

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISL) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO JK=1,KFC
    KM = D_MYMS(KMLOC)
    ISL = R_NDGNH-G_NDGLU(KM)+1
    !$ACC LOOP SEQ
    DO JGL=ISL,R_NDGNH
      IF(KM /= 0) THEN
        PSOA1(JK,JGL,KMLOC) = ZZCSTS((JK-1)+1,JGL-ISL+1,KMLOC)
      ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
        PSOA1(JK,JGL,KMLOC) = ZZCSTS((JK-1)/2+1,JGL-ISL+1,KMLOC)
      ENDIF
    ENDDO
  ENDDO
ENDDO

!$ACC END DATA

!$ACC END DATA
DEALLOCATE(ZZBS)
DEALLOCATE(ZZBA)
DEALLOCATE(ZZCSTS)
DEALLOCATE(ZZCSTA)

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEINV BARRIER')
ENDIF
CALL GSTATS(453,1)

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEINV
END MODULE LEINV_MOD
