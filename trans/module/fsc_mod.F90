MODULE FSC_MOD
CONTAINS
SUBROUTINE FSC(KF_UV,KF_SCALARS,KF_SCDERS,&
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

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_TRANS
USE TPM_DISTR
USE TPM_FIELDS
USE TPM_GEOMETRY

IMPLICIT NONE
INTEGER_M , INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS 
REAL_B , INTENT(INOUT) :: PUV(:,:)
REAL_B , INTENT(IN   ) :: PSCALAR(:,:)
REAL_B , INTENT(INOUT) :: PNSDERS(:,:)
REAL_B , INTENT(  OUT) :: PEWDERS(:,:)
REAL_B , INTENT(  OUT) :: PUVDERS(:,:)

REAL_B :: ZACHTE(D%NDGL_FS)
INTEGER_M :: IMEN(D%NDGL_FS),ISTAGTF(D%NDGL_FS)


INTEGER_M :: JGL,JLON,JF,IGLG,II,IR,JM

!     ------------------------------------------------------------------

DO JGL=1,D%NDGL_FS
  IGLG = D%NPTRLS(MYSETW)+JGL-1
  ZACHTE(JGL)  = F%RACTHE(IGLG)
  IMEN(JGL)    = G%NMEN(IGLG)
  ISTAGTF(JGL) = D%NSTAGTF(JGL)
ENDDO

#ifndef HLOMP
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JGL,JLON,JF,JM,II,IR)
#endif
DO JGL=1,D%NDGL_FS

!     ------------------------------------------------------------------

!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

  IF(KF_UV > 0) THEN
    DO JLON=ISTAGTF(JGL)+1,ISTAGTF(JGL)+2*(IMEN(JGL)+1)
      DO JF=1,2*KF_UV
        PUV(JF,JLON) = PUV(JF,JLON)*ZACHTE(JGL)
      ENDDO
    ENDDO
  ENDIF

!*      1.2      N-S DERIVATIVES

  IF(KF_SCDERS > 0)THEN
    DO JLON=ISTAGTF(JGL)+1,ISTAGTF(JGL)+2*(IMEN(JGL)+1)
      DO JF=1,KF_SCALARS
        PNSDERS(JF,JLON) = PNSDERS(JF,JLON)*ZACHTE(JGL)
      ENDDO
    ENDDO
  ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

  IF(LUVDER)THEN
    DO JM=0,IMEN(JGL)
      IR = ISTAGTF(JGL)+2*JM+1
      II = IR+1
      DO JF=1,2*KF_UV
          PUVDERS(JF,IR) = -JM*PUV(JF,II)*ZACHTE(JGL)
          PUVDERS(JF,II) =  JM*PUV(JF,IR)*ZACHTE(JGL)
      ENDDO
    ENDDO
  ENDIF

!*       2.2     SCALAR VARIABLES

  IF(KF_SCDERS > 0)THEN
    DO JM=0,IMEN(JGL)
      IR = ISTAGTF(JGL)+2*JM+1
      II = IR+1
      DO JF=1,KF_SCALARS
          PEWDERS(JF,IR) = -JM*PSCALAR(JF,II)*ZACHTE(JGL)
          PEWDERS(JF,II) =  JM*PSCALAR(JF,IR)*ZACHTE(JGL)
      ENDDO
    ENDDO
  ENDIF

ENDDO
#ifndef HLOMP
!$OMP END PARALLEL DO
#endif

!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
