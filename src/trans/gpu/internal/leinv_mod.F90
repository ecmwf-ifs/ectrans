! (C) Copyright 2000- ECMWF.
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
  
USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_FIELDS      ,ONLY : F, &
     & ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
     & IZBA,IZBS,ILDZBA,ILDZBS,ITDZBA,ITDZBS,&
     & IZCST,ILDZCA,ILDZCS,ITDZCA,ITDZCS,&
     & ZAMAX, ZSMAX
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_GEN, ONLY: NOUT
USE TPM_FLT
USE BUTTERFLY_ALG_MOD
USE CUDA_GEMM_BATCHED_MOD, ONLY: CUDA_TCGEMM_BATCHED
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM)  :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PSOA1(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PAOA1(:,:,:)

!     LOCAL
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, J1, IF, JGL,JK, J,JI, IRET
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRB) :: RRPELTMINV = 100.0_JPRB

INTEGER :: ISTAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  
  
!*       1.1      PREPARATIONS.
IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.
!$ACC DATA COPYIN (S,S%ITHRESHOLD,S%LUSEFLT,RRPELTMINV) &
!$ACC&      COPY (D,D_MYMS,R,G,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX) &
!$ACC&      PRESENT (ZAA,ZAS,IZBS,IZBA,IZCST) &
!$ACC&      PRESENT (PIA) &
!$ACC&      PRESENT (ZAMAX,ZSMAX) &
!$ACC&      COPYOUT(PSOA1,PAOA1)

!$ACC KERNELS
ISL=MAX(R_NDGNH-G_NDGLU(0)+1,1)
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(3)
DO KMLOC=1,D_NUMP
   DO JGL=ISL,R_NDGNH
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

! Initialize rescaling arrays to zero
!$ACC PARALLEL LOOP COLLAPSE(2)
DO KMLOC=1,SIZE(ZAMAX,2)
  DO JK=1,SIZE(ZAMAX,1)
    ZAMAX(JK,KMLOC) = 0.0_JPRBT
    ZSMAX(JK,KMLOC) = 0.0_JPRBT
  ENDDO
ENDDO


! 1. +++++++++++++ anti-symmetric
!$ACC PARALLEL
!$ACC LOOP
DO KMLOC=1,D_NUMP
  !$ACC LOOP SEQ
  DO J=1,(R_NSMAX+2)/2
    !$ACC LOOP PRIVATE(KM,ISKIP,ILA,IA)
    DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      IF (KM == 0) THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF

      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         ILA = (R_NSMAX-KM+2)/2
         IF (J .LE. ILA) THEN
            IA  = 1+MOD(R_NSMAX-KM+2,2)
            ZAMAX((JK-1)/ISKIP+1,KMLOC) = MAX(ZAMAX((JK-1)/ISKIP+1,KMLOC), ABS(PIA(JK,IA+1+(J-1)*2,KMLOC)))
         ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$ACC END PARALLEL

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISKIP)
DO KMLOC=1,D_NUMP
   DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      IF (KM == 0) THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF
      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         IF (ZAMAX((JK-1)/ISKIP+1,KMLOC) < 100*EPSILON(1.0_JPRB)) ZAMAX((JK-1)/ISKIP+1,KMLOC) = RRPELTMINV
      ENDIF
   ENDDO
ENDDO

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA)
DO KMLOC=1,D_NUMP
  DO J=1,(R_NSMAX+2)/2
    DO JK=1,KFC
   
      KM =  D_MYMS(KMLOC)
      IF(KM == 0)THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF
      
      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         ILA = (R_NSMAX-KM+2)/2
         IF (J .LE. ILA) THEN
            IA  = 1+MOD(R_NSMAX-KM+2,2)
            IZBA((JK-1)/ISKIP+1,J,KMLOC)=PIA(JK,IA+1+(J-1)*2,KMLOC)*RRPELTMINV/ZAMAX((JK-1)/ISKIP+1,KMLOC)
         ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO

ITHRESHOLD=S%ITHRESHOLD

!$ACC HOST_DATA USE_DEVICE(ZAA,IZBS,IZCST)
!!CALL CUDA_TCGEMM_BATCHED( &  !! the tensor core verion
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'T', &
  & ITDZCA, ILDZCA, ILDZBA, &
  & 1.0_JPRBT, &
  & IZBA, ITDZBA, ILDZBA, &
  & ZAA, LDZAA, TDZAA, &
  & 0._JPRBT, &
  & IZCST, ITDZCA, ILDZCA, &
  & D_NUMP)
!$ACC END HOST_DATA

!$ACC KERNELS 
DO KMLOC=1,D_NUMP
   DO JI=1,R_NDGNH
      DO JK=1,KFC
         KM = D_MYMS(KMLOC)
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (JI .LE. KDGLU) then
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            END IF
            
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               PAOA1(JK,ISL+JI-1,KMLOC) = IZCST((JK-1)/ISKIP+1,JI,KMLOC)*ZAMAX((JK-1)/ISKIP+1,KMLOC)/RRPELTMINV
            END IF
         END IF
      ENDDO
   ENDDO
END DO
!$ACC END KERNELS

!$ACC PARALLEL
!$ACC LOOP
DO KMLOC=1,D_NUMP
  !$ACC LOOP SEQ
  DO J=1,(R_NSMAX+3)/2
    !$ACC LOOP PRIVATE(KM,ISKIP,ILS,IS)
    DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      IF (KM == 0) THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF

      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         ILS = (R_NSMAX-KM+3)/2
         IF ((KM == 0 .AND. J .LE. ILS-1) .OR. (KM /= 0 .AND. J .LE. ILS)) THEN
            IS = 1+MOD(R_NSMAX-KM+1,2)
            ZSMAX((JK-1)/ISKIP+1,KMLOC) = MAX(ZSMAX((JK-1)/ISKIP+1,KMLOC), ABS(PIA(JK,IS+1+(J-1)*2,KMLOC)))
         END IF
      END IF
    ENDDO
  ENDDO
ENDDO
!$ACC END PARALLEL

!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISKIP)
DO KMLOC=1,D_NUMP
   DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      IF (KM == 0) THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF
      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         IF (ZSMAX((JK-1)/ISKIP+1,KMLOC) < 100*EPSILON(1.0_JPRB)) ZSMAX((JK-1)/ISKIP+1,KMLOC) = RRPELTMINV
      ENDIF
   ENDDO
ENDDO

! 2. +++++++++++++ symmetric

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILS,IS)
DO KMLOC=1,D_NUMP
  DO J=1,(R_NSMAX+3)/2
    DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      IF(KM == 0)THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF
      
      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         ILS = (R_NSMAX-KM+3)/2
         IF (J .LE. ILS) THEN
            IS  = 1+MOD(R_NSMAX-KM+1,2)
            IZBS((JK-1)/ISKIP+1,J,KMLOC)=PIA(JK,IS+1+(J-1)*2,KMLOC)*RRPELTMINV/ZSMAX((JK-1)/ISKIP+1,KMLOC)
         END IF
      END IF
    ENDDO
  ENDDO
ENDDO

!C=A*B =>
! C^T=B^T*A^T

!$ACC HOST_DATA USE_DEVICE(ZAS,IZBS,IZCST)
!!CALL CUDA_TCGEMM_BATCHED( & !! tensor-core version
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'T', &
  & ITDZCS, ILDZCS, ILDZBS, &
  & 1.0_JPRB, &
  & IZBS, ITDZBS, ILDZBS, &
  & ZAS, LDZAS, TDZAS, &
  & 0._JPRB, &
  & IZCST, ITDZCS, ILDZCS, &
  & D_NUMP)
!$ACC END HOST_DATA

!$ACC KERNELS 
DO KMLOC=1,D_NUMP
   DO JI=1,R_NDGNH
      DO JK=1,KFC
         KM = D_MYMS(KMLOC)
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (JI .LE. KDGLU) then
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            END IF
            
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               PSOA1(JK,ISL+JI-1,KMLOC) = IZCST((JK-1)/ISKIP+1,JI,KMLOC)*ZSMAX((JK-1)/ISKIP+1,KMLOC)/RRPELTMINV
            END IF
         END IF
      ENDDO
   ENDDO
END DO
!$ACC END KERNELS    

!$ACC END DATA

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEINV
END MODULE LEINV_MOD
