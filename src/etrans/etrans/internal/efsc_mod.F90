MODULE EFSC_MOD
CONTAINS
SUBROUTINE EFSC(KGL,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_TRANS       ,ONLY : LUVDER
USE TPM_DISTR       ,ONLY : D, MYSETW
!USE TPM_FIELDS
USE TPM_GEOMETRY    ,ONLY : G
USE TPMALD_GEO      ,ONLY : GALD
!

IMPLICIT NONE

INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) , INTENT(INOUT) :: PUV(:,:)
REAL(KIND=JPRB) , INTENT(IN   ) :: PSCALAR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PNSDERS(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PEWDERS(:,:)
REAL(KIND=JPRB) , INTENT(  OUT) :: PUVDERS(:,:)

INTEGER(KIND=JPIM) :: IMEN,ISTAGTF

INTEGER(KIND=JPIM) :: JF,IGLG,II,IR,JM
REAL(KIND=JPRB) :: ZIM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EFSC_MOD:EFSC',0,ZHOOK_HANDLE)
IGLG    = D%NPTRLS(MYSETW)+KGL-1
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)

!     ------------------------------------------------------------------

!*           EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  DO JM=0,IMEN
    ZIM=REAL(JM,JPRB)*GALD%EXWN
    IR = ISTAGTF+2*JM+1
    II = IR+1
! use unroll to provoke vectorization of outer loop
!cdir unroll=4
    DO JF=1,2*KF_UV
      PUVDERS(JF,IR) = -PUV(JF,II)*ZIM
      PUVDERS(JF,II) =  PUV(JF,IR)*ZIM
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
  DO JM=0,IMEN
    ZIM=REAL(JM,JPRB)*GALD%EXWN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    DO JF=1,KF_SCALARS
      PEWDERS(JF,IR) = -PSCALAR(JF,II)*ZIM
      PEWDERS(JF,II) =  PSCALAR(JF,IR)*ZIM
    ENDDO
  ENDDO
ENDIF
IF (LHOOK) CALL DR_HOOK('EFSC_MOD:EFSC',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EFSC
END MODULE EFSC_MOD
