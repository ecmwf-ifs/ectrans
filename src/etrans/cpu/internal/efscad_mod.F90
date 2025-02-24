! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EFSCAD_MOD
CONTAINS
SUBROUTINE EFSCAD(KGL,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

!**** *EFSCAD - Division by a*cos(theta), east-west derivatives - adjoint

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL EFSCAD(..)
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

IMPLICIT NONE

INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) , INTENT(INOUT) :: PUV(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PSCALAR(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PNSDERS(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PEWDERS(:,:)
REAL(KIND=JPRB) , INTENT(INOUT) :: PUVDERS(:,:)

INTEGER(KIND=JPIM) :: IMEN,ISTAGTF

INTEGER(KIND=JPIM) :: JF,IGLG,II,IR,JM

REAL(KIND=JPRB) :: ZIM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EFSCAD_MOD:EFSCAD',0,ZHOOK_HANDLE)
IGLG = D%NPTRLS(MYSETW)+KGL-1
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  DO JM=0,IMEN

    ZIM=REAL(JM,JPRB)*GALD%EXWN

    IR = ISTAGTF+2*JM+1
    II = IR+1
    DO JF=1,2*KF_UV

      PUV(JF,II) = PUV(JF,II) - ZIM*PUVDERS(JF,IR)
      PUV(JF,IR) = PUV(JF,IR) + ZIM*PUVDERS(JF,II)

      PUVDERS(JF,IR) = 0.0_JPRB
      PUVDERS(JF,II) = 0.0_JPRB
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

      PSCALAR(JF,II) = PSCALAR(JF,II) - ZIM* PEWDERS(JF,IR)
      PSCALAR(JF,IR) = PSCALAR(JF,IR) + ZIM* PEWDERS(JF,II)

      PEWDERS(JF,IR) = 0.0_JPRB
      PEWDERS(JF,II) = 0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF (LHOOK) CALL DR_HOOK('EFSCAD_MOD:EFSCAD',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE EFSCAD
END MODULE EFSCAD_MOD
