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
CONTAINS
SUBROUTINE LEDIR(KM,KMLOC,KFC,KIFC,KSL,KDGLU,KLED2,PAIA,PSIA,POA1,PW)

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
INTEGER(KIND=JPIM), INTENT(IN)  :: KLED2

REAL(KIND=JPRD),    INTENT(IN)  :: PW(KDGLU+KSL-1)
REAL(KIND=JPRB),    INTENT(IN)  :: PSIA(:,:),   PAIA(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: POA1(:,:)

!     LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IFLD, J, JK, I1, I2, I3, I4
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRB)    :: ZB(KDGLU,KIFC), ZCA((R%NTMAX-KM+2)/2,KIFC), ZCS((R%NTMAX-KM+3)/2,KIFC)
REAL(KIND=JPRD), allocatable :: ZB_D(:,:), ZCA_D(:,:), ZCS_D(:,:),ZRPNMA(:,:), ZRPNMS(:,:)
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
ISL = KSL

IF(KM == 0)THEN
  ISKIP = 2
ELSE
  ISKIP = 1
ENDIF

IF (KIFC > 0 .AND. KDGLU > 0 ) THEN

  ITHRESHOLD=S%ITHRESHOLD
 
!*       1. ANTISYMMETRIC PART.

  IFLD=0
  DO JK=1,KFC,ISKIP
    IFLD=IFLD+1
    DO J=1,KDGLU
      ZB(J,IFLD)=PAIA(JK,ISL+J-1)*PW(ISL+J-1)
    ENDDO
  ENDDO
  
  IF(ILA <= ITHRESHOLD .OR. .NOT.S%LUSEFLT) THEN

    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1',0,ZHOOK_HANDLE)
    IF (LLDOUBLE) THEN
       CALL GEMM('T','N',ILA,KIFC,KDGLU,1.0_JPRD,S%FA(KMLOC)%RPNMA,KDGLU,&
            &ZB,KDGLU,0._JPRD,ZCA,ILA)
    ELSE
       IF(KM>=1)THEN ! DGEM for the mean to improve mass conservation
          CALL GEMM('T','N',ILA,KIFC,KDGLU,1.0_JPRM,S%FA(KMLOC)%RPNMA,KDGLU,&
               &ZB,KDGLU,0._JPRM,ZCA,ILA)
       ELSE
          I1 = size(S%FA(KMLOC)%RPNMA(:,1))
          I2 = size(S%FA(KMLOC)%RPNMA(1,:))
          ALLOCATE(ZRPNMA(I1,I2))
          ALLOCATE(ZB_D(KDGLU,KIFC))
          ALLOCATE(ZCA_D((R%NTMAX-KM+2)/2,KIFC))
          IFLD=0
          DO JK=1,KFC,ISKIP
             IFLD=IFLD+1
             DO J=1,KDGLU
                ZB_D(J,IFLD)=ZB(J,IFLD)
             ENDDO
          ENDDO
          DO I3=1,I1
             DO I4=1,I2
                ZRPNMA(I3,I4) = S%FA(KMLOC)%RPNMA(I3,I4)
             END DO
          END DO
          CALL GEMM('T','N',ILA,KIFC,KDGLU,1.0_JPRD,ZRPNMA,KDGLU,&
               &ZB_D,KDGLU,0._JPRD,ZCA_D,ILA)
          IFLD=0
          DO JK=1,KFC,ISKIP
             IFLD=IFLD+1
             DO J=1,ILA
                ZCA(J,IFLD) = ZCA_D(J,IFLD)
             ENDDO
          ENDDO
          DEALLOCATE(ZRPNMA)
          DEALLOCATE(ZB_D)
          DEALLOCATE(ZCA_D)
       END IF
    ENDIF
    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_1',1,ZHOOK_HANDLE)

  ELSE
     IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_1',0,ZHOOK_HANDLE)
     CALL MULT_BUTM('T',S%FA(KMLOC)%YBUT_STRUCT_A,KIFC,ZB,ZCA,KM)
     IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_1',1,ZHOOK_HANDLE)
  ENDIF

  IFLD=0
  DO JK=1,KFC,ISKIP
    IFLD=IFLD+1
    DO J=1,ILA
      POA1(IA+(J-1)*2,JK) = ZCA(J,IFLD)
    ENDDO
  ENDDO
  
!*       1.3      SYMMETRIC PART.

  
  IFLD=0
  DO JK=1,KFC,ISKIP
    IFLD=IFLD+1
    DO J=1,KDGLU
      ZB(J,IFLD)=PSIA(JK,ISL+J-1)*REAL(PW(ISL+J-1),JPRB)
    ENDDO
  ENDDO
  
  IF(ILS <= ITHRESHOLD .OR. .NOT.S%LUSEFLT) THEN

    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2',0,ZHOOK_HANDLE)
    IF (LLDOUBLE) THEN
       CALL GEMM('T','N',ILS,KIFC,KDGLU,1.0_JPRD,S%FA(KMLOC)%RPNMS,KDGLU,&
            &ZB,KDGLU,0._JPRD,ZCS,ILS)
    ELSE
       IF(KM>=1)THEN ! DGEM for the mean to improve mass conservation
          CALL GEMM('T','N',ILS,KIFC,KDGLU,1.0_JPRM,S%FA(KMLOC)%RPNMS,KDGLU,&
               &ZB,KDGLU,0._JPRM,ZCS,ILS)
       ELSE
          I1 = size(S%FA(KMLOC)%RPNMS(:,1))
          I2 = size(S%FA(KMLOC)%RPNMS(1,:))
          ALLOCATE(ZRPNMS(I1,I2))
          ALLOCATE(ZB_D(KDGLU,KIFC))
          ALLOCATE(ZCS_D((R%NTMAX-KM+3)/2,KIFC))          
          IFLD=0
          DO JK=1,KFC,ISKIP
             IFLD=IFLD+1
             DO J=1,KDGLU
                ZB_D(J,IFLD)=PSIA(JK,ISL+J-1)*REAL(PW(ISL+J-1),JPRB)
             ENDDO
          ENDDO
          DO I3=1,I1
             DO I4=1,I2
                ZRPNMS(I3,I4) = S%FA(KMLOC)%RPNMS(I3,I4)
             END DO
          END DO
          CALL GEMM('T','N',ILS,KIFC,KDGLU,1.0_JPRD,ZRPNMS,KDGLU,&
               &ZB_D,KDGLU,0._JPRD,ZCS_D,ILS)
          IFLD=0
          DO JK=1,KFC,ISKIP
             IFLD=IFLD+1
             DO J=1,ILS
                ZCS(J,IFLD) = ZCS_D(J,IFLD)
             ENDDO
          ENDDO
          DEALLOCATE(ZRPNMS)
          DEALLOCATE(ZB_D)
          DEALLOCATE(ZCS_D)
       END IF
    ENDIF
    IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'GEMM_2',1,ZHOOK_HANDLE)
    
  ELSE
     IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_2',0,ZHOOK_HANDLE)
     CALL MULT_BUTM('T',S%FA(KMLOC)%YBUT_STRUCT_S,KIFC,ZB,ZCS,KM)
     IF (LHOOK) CALL DR_HOOK('LEDIR_'//CLX//'BUTM_2',1,ZHOOK_HANDLE)
  ENDIF

  IFLD=0
  DO JK=1,KFC,ISKIP
    IFLD=IFLD+1
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
