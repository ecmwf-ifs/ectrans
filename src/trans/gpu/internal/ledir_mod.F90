! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEDIR_MOD
CONTAINS
SUBROUTINE LEDIR(KF_FS,KLED2,PAIA,PSIA,POA1)

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

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX,R_NTMAX
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_FIELDS      ,ONLY : F, &
     & ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
     & DZBST,DLDZBA,DLDZBS,DTDZBA,DTDZBS,&
     & DZCST,DZCAT,DLDZCA,DLDZCS,DTDZCA,DTDZCS,&
     & ZAMAX, ZSMAX
USE TPM_DISTR
USE TPM_GEN, ONLY: NOUT
USE TPM_FLT
USE BUTTERFLY_ALG_MOD
USE CUDA_GEMM_BATCHED_MOD, ONLY: CUDA_TCGEMM_BATCHED, CUDA_GEMM_BATCHED
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM)  :: KFC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM) :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN)  :: KLED2

REAL(KIND=JPRBT),    INTENT(IN)  :: PSIA(:,:,:),   PAIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: POA1(:,:,:)

!     LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IF, J, JK, IRET
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRB) :: RRPELTMDIR = 100.0_JPRB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER :: ISTAT

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)

KFC = 2*KF_FS
KIFC = KFC

!$ACC DATA &
!$ACC& COPY(F,F%RW) &
!$ACC& COPY(D,D_NUMP,D_MYMS,R,R_NDGNH,G,G_NDGLU,R_NSMAX,R_NTMAX) &
!$ACC& PRESENT(PSIA,PAIA) &
!$ACC& PRESENT(ZAA,ZAS,DZBST,DZCST,DZCAT) &
!$ACC& PRESENT (ZAMAX,ZSMAX) &
!$ACC& PRESENT(POA1)

! Initialize rescaling arrays to zero
!$ACC PARALLEL LOOP COLLAPSE(2)
DO KMLOC=1,SIZE(ZAMAX,2)
  DO JK=1,SIZE(ZAMAX,1)
    ZAMAX(JK,KMLOC) = 0.0_JPRBT
    ZSMAX(JK,KMLOC) = 0.0_JPRBT
  ENDDO
ENDDO

! anti-symmetric
!$ACC PARALLEL
!$ACC LOOP
DO KMLOC=1,D_NUMP
  !$ACC LOOP SEQ
  DO J=1,R_NDGNH
    !$ACC LOOP PRIVATE(KM,KDGLU,ISKIP)
    DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
      IF (J .LE. KDGLU) THEN
        ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
        IF (KM == 0) THEN
           ISKIP = 2
        ELSE
           ISKIP = 1
        ENDIF

        IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
          ZAMAX((JK-1)/ISKIP+1,KMLOC) = MAX(ZAMAX((JK-1)/ISKIP+1,KMLOC), ABS(PAIA(JK,ISL+J-1,KMLOC)))
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
         IF (ZAMAX((JK-1)/ISKIP+1,KMLOC) < 100*EPSILON(1.0_JPRB)) ZAMAX((JK-1)/ISKIP+1,KMLOC) = RRPELTMDIR
      ENDIF
   ENDDO
ENDDO

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP)
DO KMLOC=1,D_NUMP
   DO J=1,R_NDGNH   
      DO JK=1,KFC
         
         KM = D_MYMS(KMLOC)   
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (J .LE. KDGLU) THEN
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            ENDIF
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               DZBST((JK-1)/ISKIP+1,J,KMLOC)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)*RRPELTMDIR/ZAMAX((JK-1)/ISKIP+1,KMLOC)
            END IF
         END IF
      ENDDO
   ENDDO
END DO


!$ACC HOST_DATA USE_DEVICE(ZAA,DZBST,DZCAT)
!!CALL CUDA_TCGEMM_BATCHED( &  !! tensor-core version
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'N', &
  & DTDZBA, TDZAA, DLDZBA, &
  & 1.0_JPRBT, &
  & DZBST, DTDZBA, DLDZBA, &
  & ZAA, LDZAA, TDZAA, &
  & 0._JPRBT, &
  & DZCAT, DTDZCA, DLDZCA, &
  & D_NUMP)
!$ACC END HOST_DATA

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA)
DO KMLOC=1,D_NUMP
   DO J=1,(R_NTMAX+2)/2
      DO JK=1,KFC

         KM = D_MYMS(KMLOC)
         IF(KM == 0)THEN
            ISKIP = 2
         ELSE
            ISKIP = 1
         ENDIF
         
         IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            ILA = (R_NTMAX-KM+2)/2
            IA  = 1+MOD(R_NTMAX-KM+2,2)
            IF (J .LE. ILA) THEN
               POA1(JK,IA+(J-1)*2,KMLOC) = DZCAT((JK-1)/ISKIP+1,J,KMLOC)*ZAMAX((JK-1)/ISKIP+1,KMLOC)/RRPELTMDIR
            END IF
         END IF
      ENDDO
   ENDDO
ENDDO

! symmetric
!$ACC PARALLEL
!$ACC LOOP
DO KMLOC=1,D_NUMP
  !$ACC LOOP SEQ
  DO J=1,R_NDGNH
    !$ACC LOOP PRIVATE(KM,KDGLU,ISKIP)
    DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
      IF (J .LE. KDGLU) THEN
        IF (KM == 0) THEN
           ISKIP = 2
        ELSE
           ISKIP = 1
        ENDIF

        IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
          ZSMAX((JK-1)/ISKIP+1,KMLOC) = MAX(ZSMAX((JK-1)/ISKIP+1,KMLOC), ABS(PSIA(JK,ISL+J-1,KMLOC)))
        END IF
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
         IF (ZSMAX((JK-1)/ISKIP+1,KMLOC) < 100*EPSILON(1.0_JPRB)) ZSMAX((JK-1)/ISKIP+1,KMLOC) = RRPELTMDIR
      ENDIF
   ENDDO
ENDDO

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP)
DO KMLOC=1,D_NUMP
   DO J=1,R_NDGNH   
      DO JK=1,KFC
         KM = D_MYMS(KMLOC)   
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (J .LE. KDGLU) THEN
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            ENDIF
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               DZBST((JK-1)/ISKIP+1,J,KMLOC)=PSIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)*RRPELTMDIR/ZSMAX((JK-1)/ISKIP+1,KMLOC)
            END IF
         END IF
      ENDDO
   ENDDO
END DO

!$ACC HOST_DATA USE_DEVICE(ZAS,DZBST,DZCST)
!!CALL CUDA_TCGEMM_BATCHED( & !! tensor-core version
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'N', &
  & DTDZBS, TDZAS, DLDZBS, &
  & 1.0_JPRBT, &
  & DZBST, DTDZBS, DLDZBS, &
  & ZAS, LDZAS, TDZAS, &
  & 0._JPRBT, &
  & DZCST, DTDZCS, DLDZCS, &
  & D_NUMP)
!$ACC END HOST_DATA
!    CHARACTER,                         INTENT(IN)  :: TRANSA
!    CHARACTER,                         INTENT(IN)  :: TRANSB
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: M
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: N
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: K
!    REAL(KIND=JPRD),                   INTENT(IN)  :: ALPHA
!    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(IN)  :: AARRAY
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDA
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEA
!    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(IN)  :: BARRAY
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDB
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEB
!    REAL(KIND=JPRD),                   INTENT(IN)  :: BETA
!    REAL(KIND=JPRD), DIMENSION(:,:,:), INTENT(OUT) :: CARRAY
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: LDC
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: STRIDEC
!    INTEGER(KIND=JPIM),                INTENT(IN)  :: BATCHCOUNT

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA,ILS,IS)
DO KMLOC=1,D_NUMP
   DO J=1,(R_NTMAX+3)/2
      DO JK=1,KFC

         KM = D_MYMS(KMLOC)
         IF(KM == 0)THEN
            ISKIP = 2
         ELSE
            ISKIP = 1
         ENDIF
         
         IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            ILS = (R_NTMAX-KM+3)/2
            IF (J .LE. ILS) THEN
               IS  = 1+MOD(R_NTMAX-KM+1,2)
               POA1(JK,IS+(J-1)*2,KMLOC) = DZCST((JK-1)/ISKIP+1,J,KMLOC)*ZSMAX((JK-1)/ISKIP+1,KMLOC)/RRPELTMDIR
            END IF
         END IF
      ENDDO
   ENDDO
ENDDO

!$ACC END DATA

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
