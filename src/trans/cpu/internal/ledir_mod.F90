! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEDIR_MOD

PRIVATE
PUBLIC LEDIR

CONTAINS
SUBROUTINE LEDIR(KM,KMLOC,KFC,KIFC,KSL,KDGLU,PAIA,PSIA,POA1,PW)

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
!      P. Dueben : Dec 2019 Improvements for mass conservation in single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPRM, JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_FLT         ,ONLY : S
USE BUTTERFLY_ALG_MOD, ONLY : MULT_BUTM
USE ECTRANS_BLAS_MOD, ONLY : GEMM

IMPLICIT NONE

!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KIFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KSL
INTEGER(KIND=JPIM), INTENT(IN)  :: KDGLU

REAL(KIND=JPRD),    INTENT(IN)  :: PW(KDGLU+KSL-1)
REAL(KIND=JPRB),    INTENT(IN)  :: PSIA(:,:),   PAIA(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: POA1(:,:)

!     LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, IFLD, J, JK
REAL(KIND=JPRB)    :: ZB(KDGLU,KIFC), ZCA((R%NTMAX-KM+2)/2,KIFC), ZCS((R%NTMAX-KM+3)/2,KIFC)
LOGICAL, PARAMETER :: LLDOUBLE = (JPRB == JPRD)
CHARACTER(LEN=1) :: CLX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

CLX = 'S'
IF (LLDOUBLE) CLX = 'D'

IA  = 1+MOD(R%NTMAX-KM+2,2)
IS  = 1+MOD(R%NTMAX-KM+1,2)
ILA = (R%NTMAX-KM+2)/2
ILS = (R%NTMAX-KM+3)/2

ISKIP = MERGE(2, 1, KM == 0)

IF (KIFC > 0 .AND. KDGLU > 0 ) THEN
  !*       1. ANTISYMMETRIC PART.
  IF (ILA <= S%ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    IF (LLDOUBLE) THEN
      CALL PACK_FOR_GEMM(ISKIP, KFC, KIFC, KDGLU, KSL, PAIA, PW, ZB)

      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 0, ZHOOK_HANDLE)
      CALL GEMM('T', 'N', ILA, KIFC, KDGLU, 1.0_JPRD, S%FA(KMLOC)%RPNMA, KDGLU, ZB, KDGLU, &
        &       0._JPRD, ZCA, ILA)
      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 1, ZHOOK_HANDLE)
    ELSE
      IF (KM >= 1) THEN
        CALL PACK_FOR_GEMM(ISKIP, KFC, KIFC, KDGLU, KSL, PAIA, PW, ZB)

        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 0, ZHOOK_HANDLE)
        CALL GEMM('T', 'N', ILA, KIFC, KDGLU, 1.0_JPRM, S%FA(KMLOC)%RPNMA, KDGLU, ZB, KDGLU, &
              &   0._JPRM, ZCA, ILA)
        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 1, ZHOOK_HANDLE)
      ELSE
        BLOCK
          REAL(KIND=JPRD), ALLOCATABLE :: ZB_D(:,:), ZCA_D(:,:)

          ALLOCATE(ZB_D(KDGLU,KIFC), ZCA_D(ILA,KIFC))

          CALL PACK_FOR_GEMM_JPRD(ISKIP, KFC, KIFC, KDGLU, KSL, PAIA, PW, ZB_D)

          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 0, ZHOOK_HANDLE)
          CALL GEMM('T', 'N', ILA, KIFC, KDGLU, 1.0_JPRD, S%RPNMA_DGEMM, KDGLU, ZB_D, KDGLU, &
              &       0._JPRD, ZCA_D, ILA)
          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 1, ZHOOK_HANDLE)

          ZCA(:,:) = ZCA_D(:,:)

          DEALLOCATE(ZB_D, ZCA_D)
        END BLOCK
      END IF
    ENDIF
  ELSE

    CALL PACK_FOR_GEMM(ISKIP, KFC, KIFC, KDGLU, KSL, PAIA, PW, ZB)

    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_1', 0, ZHOOK_HANDLE)
    CALL MULT_BUTM('T', S%FA(KMLOC)%YBUT_STRUCT_A, KIFC, ZB, ZCA, KM)
    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_1', 1, ZHOOK_HANDLE)
  ENDIF

  !*       1.3      SYMMETRIC PART.
  IF (ILS <= S%ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    IF (LLDOUBLE) THEN
      CALL PACK_FOR_GEMM(ISKIP, KFC, KIFC, KDGLU, KSL, PSIA, PW, ZB)

      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 0, ZHOOK_HANDLE)
      CALL GEMM('T', 'N', ILS, KIFC, KDGLU, 1.0_JPRD, S%FA(KMLOC)%RPNMS, KDGLU, ZB, KDGLU, &
        &       0._JPRD, ZCS, ILS)
      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 1, ZHOOK_HANDLE)
    ELSE
      IF (KM >= 1) THEN
        CALL PACK_FOR_GEMM(ISKIP, KFC, KIFC, KDGLU, KSL, PSIA, PW, ZB)

        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 0, ZHOOK_HANDLE)
        CALL GEMM('T', 'N', ILS, KIFC, KDGLU, 1.0_JPRM, S%FA(KMLOC)%RPNMS, KDGLU, ZB, KDGLU, &
          &       0._JPRM, ZCS, ILS)
        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 1, ZHOOK_HANDLE)
      ELSE
        BLOCK
          REAL(KIND=JPRD), ALLOCATABLE :: ZB_D(:,:), ZCS_D(:,:)

          ALLOCATE(ZB_D(KDGLU,KIFC), ZCS_D(ILS,KIFC))

          CALL PACK_FOR_GEMM_JPRD(ISKIP, KFC, KIFC, KDGLU, KSL, PSIA, PW, ZB_D)

          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 0, ZHOOK_HANDLE)
          CALL GEMM('T', 'N', ILS, KIFC, KDGLU, 1.0_JPRD, S%RPNMS_DGEMM, KDGLU, ZB_D, KDGLU, &
            &       0._JPRD, ZCS_D, ILS)
          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 1, ZHOOK_HANDLE)

          ZCS(:,:) = ZCS_D(:,:)

          DEALLOCATE(ZB_D, ZCS_D)
        END BLOCK
      END IF
    ENDIF
  ELSE
    CALL PACK_FOR_GEMM(ISKIP, KFC, KIFC, KDGLU, KSL, PSIA, PW, ZB)

    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_2', 0, ZHOOK_HANDLE)
    CALL MULT_BUTM('T', S%FA(KMLOC)%YBUT_STRUCT_S, KIFC, ZB, ZCS, KM)
    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_2', 1, ZHOOK_HANDLE)
  ENDIF

  IFLD=0
  DO JK=1,KFC,ISKIP
    IFLD=IFLD+1
    DO J=1,ILA
      POA1(IA+(J-1)*2,JK) = ZCA(J,IFLD)
    ENDDO
    DO J=1,ILS
      POA1(IS+(J-1)*2,JK) = ZCS(J,IFLD)
    ENDDO
  ENDDO
ELSE
  ! This zonal wavenumber KM has no computation to be done (G%NDGLU(KM) = 0)
  ! This is usually because the wavenumber cannot be represented on the given grid, so we should
  ! zero POA1
  POA1(:,:) = 0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LEDIR

SUBROUTINE PACK_FOR_GEMM(KSKIP, KFC, KIFC, KDGLU, KSL, PIN, PW, POUT)

USE PARKIND1, ONLY: JPIM, JPRB, JPRD

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KSKIP
INTEGER(KIND=JPIM), INTENT(IN) :: KFC
INTEGER(KIND=JPIM), INTENT(IN) :: KIFC
INTEGER(KIND=JPIM), INTENT(IN) :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN) :: KSL
REAL(KIND=JPRB), INTENT(IN) :: PIN(:, :)
REAL(KIND=JPRD), INTENT(IN) :: PW(:)
REAL(KIND=JPRB), INTENT(OUT) :: POUT(KDGLU, KIFC)

INTEGER(KIND=JPIM) :: IFLD, JK, J

IFLD = 0
DO JK = 1, KFC, KSKIP
  IFLD = IFLD + 1
  DO J = 1, KDGLU
    POUT(J, IFLD) = PIN(JK, J) * REAL(PW(KSL + J - 1), JPRB)
  ENDDO
ENDDO

END SUBROUTINE PACK_FOR_GEMM

SUBROUTINE PACK_FOR_GEMM_JPRD(KSKIP, KFC, KIFC, KDGLU, KSL, PIN, PW, POUT)

USE PARKIND1, ONLY: JPIM, JPRB, JPRD

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KSKIP
INTEGER(KIND=JPIM), INTENT(IN) :: KFC
INTEGER(KIND=JPIM), INTENT(IN) :: KIFC
INTEGER(KIND=JPIM), INTENT(IN) :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN) :: KSL
REAL(KIND=JPRB), INTENT(IN) :: PIN(:, :)
REAL(KIND=JPRD), INTENT(IN) :: PW(:)
REAL(KIND=JPRD), INTENT(OUT) :: POUT(KDGLU, KIFC)

INTEGER(KIND=JPIM) :: IFLD, JK, J

IFLD = 0
DO JK = 1, KFC, KSKIP
  IFLD = IFLD + 1
  DO J = 1, KDGLU
    POUT(J, IFLD) = REAL(PIN(JK, J), JPRD) * PW(KSL + J - 1)
  ENDDO
ENDDO

END SUBROUTINE PACK_FOR_GEMM_JPRD

END MODULE LEDIR_MOD
