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
REAL(KIND=JPRB),    INTENT(INOUT)  :: PSIA(:,:),   PAIA(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: POA1(:,:)

!     LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, IFLD, J, JK
REAL(KIND=JPRB)    :: ZCA((R%NTMAX-KM+2)/2,KIFC), ZCS((R%NTMAX-KM+3)/2,KIFC)
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
  DO IFLD = 1, KIFC
    PAIA(:, IFLD) = REAL(REAL(PAIA(:, IFLD), JPRD) * PW(KSL:KSL + KDGLU - 1), JPRB)
  ENDDO

  IF (ILA <= S%ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    IF (LLDOUBLE) THEN
      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 0, ZHOOK_HANDLE)
      CALL GEMM('T', 'N', ILA, KIFC, KDGLU, 1.0_JPRD, S%FA(KMLOC)%RPNMA, KDGLU, PAIA, KDGLU, &
        &       0._JPRD, ZCA, ILA)
      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 1, ZHOOK_HANDLE)
    ELSE
      IF (KM >= 1) THEN
        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 0, ZHOOK_HANDLE)
        CALL GEMM('T', 'N', ILA, KIFC, KDGLU, 1.0_JPRM, S%FA(KMLOC)%RPNMA, KDGLU, PAIA, KDGLU, &
              &   0._JPRM, ZCA, ILA)
        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 1, ZHOOK_HANDLE)
      ELSE
        BLOCK
          REAL(KIND=JPRD), ALLOCATABLE :: ZAIA_D(:,:), ZCA_D(:,:)

          ALLOCATE(ZAIA_D(KDGLU,KIFC), ZCA_D(ILA,KIFC))
          ZAIA_D(:,:) = REAL(PAIA(:,:), JPRD)

          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 0, ZHOOK_HANDLE)
          CALL GEMM('T', 'N', ILA, KIFC, KDGLU, 1.0_JPRD, S%RPNMA_DGEMM, KDGLU, ZAIA_D, KDGLU, &
              &       0._JPRD, ZCA_D, ILA)
          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1', 1, ZHOOK_HANDLE)

          ZCA(:,:) = ZCA_D(:,:)

          DEALLOCATE(ZAIA_D, ZCA_D)
        END BLOCK
      END IF
    ENDIF
  ELSE
    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_1', 0, ZHOOK_HANDLE)
    CALL MULT_BUTM('T', S%FA(KMLOC)%YBUT_STRUCT_A, KIFC, PAIA, ZCA, KM)
    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_1', 1, ZHOOK_HANDLE)
  ENDIF

  !*       1.3      SYMMETRIC PART.
  DO IFLD = 1, KIFC
    PSIA(:, IFLD) = REAL(REAL(PSIA(:, IFLD), JPRD) * PW(KSL:KSL + KDGLU - 1), JPRB)
  ENDDO

  IF (ILS <= S%ITHRESHOLD .OR. .NOT. S%LUSEFLT) THEN
    IF (LLDOUBLE) THEN
      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 0, ZHOOK_HANDLE)
      CALL GEMM('T', 'N', ILS, KIFC, KDGLU, 1.0_JPRD, S%FA(KMLOC)%RPNMS, KDGLU, PSIA, KDGLU, &
        &       0._JPRD, ZCS, ILS)
      IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 1, ZHOOK_HANDLE)
    ELSE
      IF (KM >= 1) THEN
        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 0, ZHOOK_HANDLE)
        CALL GEMM('T', 'N', ILS, KIFC, KDGLU, 1.0_JPRM, S%FA(KMLOC)%RPNMS, KDGLU, PSIA, KDGLU, &
          &       0._JPRM, ZCS, ILS)
        IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 1, ZHOOK_HANDLE)
      ELSE
        BLOCK
          REAL(KIND=JPRD), ALLOCATABLE :: ZSIA_D(:,:), ZCS_D(:,:)

          ALLOCATE(ZSIA_D(KDGLU,KIFC), ZCS_D(ILS,KIFC))
          ZSIA_D(:,:) = REAL(PSIA(:,:), JPRD)

          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 0, ZHOOK_HANDLE)
          CALL GEMM('T', 'N', ILS, KIFC, KDGLU, 1.0_JPRD, S%RPNMS_DGEMM, KDGLU, ZSIA_D, KDGLU, &
            &       0._JPRD, ZCS_D, ILS)
          IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2', 1, ZHOOK_HANDLE)

          ZCS(:,:) = ZCS_D(:,:)

          DEALLOCATE(ZSIA_D, ZCS_D)
        END BLOCK
      END IF
    ENDIF
  ELSE
    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_2', 0, ZHOOK_HANDLE)
    CALL MULT_BUTM('T', S%FA(KMLOC)%YBUT_STRUCT_S, KIFC, PSIA, ZCS, KM)
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

END MODULE LEDIR_MOD
