! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSC_VIEW_MOD
CONTAINS
SUBROUTINE FSC_VIEW(KGL,&
 & YDPUV,YDPSCALAR,&
 & YDPUV_EW,YDPSCALAR_NS,YDPSCALAR_EW)

!**** *FSC_VIEW - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC_VIEW(..)
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
!        Original : 00-03-03 (From SC2FSC_VIEW)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD ,ONLY : SPEC_VIEW
!

IMPLICIT NONE


INTEGER(KIND=JPIM) , INTENT(IN) :: KGL
TYPE(SPEC_VIEW), INTENT(IN) :: YDPUV(:),YDPSCALAR(:)
TYPE(SPEC_VIEW), INTENT(IN) :: YDPUV_EW(:),YDPSCALAR_NS(:), YDPSCALAR_EW(:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM
INTEGER(KIND=JPIM) :: KF_UV,KF_SCALARS,KF_SCDERS

KF_UV = SIZE(YDPUV)/2
KF_SCALARS = SIZE(YDPSCALAR)
KF_SCDERS = SIZE(YDPSCALAR_NS)

!------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = REAL(F%RACTHE(IGLG),JPRB)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2 = REAL(F%RACTHE(IGLG),JPRB)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = REAL(F%RACTHE2(IGLG),JPRB)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    DO JF=1,KF_SCALARS
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        
        ! calculate amplitude and add phase shift then reconstruct A,B
        ZAMP = SQRT(YDPSCALAR(JF)%P(IR)**2 + YDPSCALAR(JF)%P(II)**2)
        ZPHASE = ATAN2(YDPSCALAR(JF)%P(II),YDPSCALAR(JF)%P(IR)) + REAL(JM,JPRB)*ZSHIFT
        
        YDPSCALAR(JF)%P(IR) = ZAMP*COS(ZPHASE)
        YDPSCALAR(JF)%P(II) = ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
    IF(KF_SCDERS > 0)THEN
      DO JF=1,KF_SCALARS
        DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(YDPSCALAR_NS(JF)%P(IR)**2 + YDPSCALAR_NS(JF)%P(II)**2)
          ZPHASE = ATAN2(YDPSCALAR_NS(JF)%P(II),YDPSCALAR_NS(JF)%P(IR)) + REAL(JM,JPRB)*ZSHIFT
          YDPSCALAR_NS(JF)%P(IR) = ZAMP*COS(ZPHASE)
          YDPSCALAR_NS(JF)%P(II) = ZAMP*SIN(ZPHASE)
        ENDDO
      ENDDO
    ENDIF
    DO JF=1,2*KF_UV
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        ! calculate amplitude and phase shift and reconstruct A,B
        ZAMP = SQRT(YDPUV(JF)%P(IR)**2 + YDPUV(JF)%P(II)**2)
        ZPHASE = ATAN2(YDPUV(JF)%P(II),YDPUV(JF)%P(IR)) + REAL(JM,JPRB)*ZSHIFT
        YDPUV(JF)%P(IR) =  ZAMP*COS(ZPHASE)
        YDPUV(JF)%P(II) =  ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------
  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,2*KF_UV
      YDPUV(JF)%P(JLON) = YDPUV(JF)%P(JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,KF_SCALARS
      YDPSCALAR_NS(JF)%P(JLON) = YDPSCALAR_NS(JF)%P(JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,2*KF_UV
      YDPUV_EW(JF)%P(IR) = -YDPUV(JF)%P(II)*ZMUL
      YDPUV_EW(JF)%P(II) =  YDPUV(JF)%P(IR)*ZMUL
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,KF_SCALARS
      YDPSCALAR_EW(JF)%P(IR) = -YDPSCALAR(JF)%P(II)*ZMUL
      YDPSCALAR_EW(JF)%P(II) =  YDPSCALAR(JF)%P(IR)*ZMUL
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSC_VIEW
END MODULE FSC_VIEW_MOD
