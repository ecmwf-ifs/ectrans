! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSC_MOD
CONTAINS
SUBROUTINE FSC(KF_UV,KF_SCALARS,KF_SCDERS,&
 & KST_UV,KST_SC,KST_NSDERS,KST_EWDERS,KST_UVDERS)

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

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT

USE TPM_TRANS       ,ONLY : LUVDER, LATLON, ZGTF
USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
use tpm_gen, only: nout
!

IMPLICIT NONE
INTEGER(KIND=JPIM) :: KGL
INTEGER(KIND=JPIM) , INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS
INTEGER(KIND=JPIM) , INTENT(IN) :: KST_UV, KST_SC, KST_NSDERS, KST_EWDERS, KST_UVDERS

REAL(KIND=JPRBT) , POINTER :: PUV(:,:)
REAL(KIND=JPRBT) , POINTER :: PSCALAR(:,:)
REAL(KIND=JPRBT) , POINTER :: PNSDERS(:,:)
REAL(KIND=JPRBT) , POINTER :: PEWDERS(:,:)
REAL(KIND=JPRBT) , POINTER :: PUVDERS(:,:)

REAL(KIND=JPRBT) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRBT) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
!DEBUGGING:
integer :: i,J,maxi,maxj
real :: maxv

!     ------------------------------------------------------------------

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

!write(301,*) ' nums ', KST_UV, KST_SC, KF_UV, KST_nsders, KST_ewders, KF_SCDERS, KST_uvders, D%NDGL_FS
IF( KF_UV > 0 ) THEN
  PUV     => ZGTF(2*KST_UV-1:2*(KST_UV+2*KF_UV-1),:)
ENDIF
PSCALAR => ZGTF(2*KST_SC-1:2*(KST_SC+KF_SCALARS-1),:)
IF( KF_SCDERS > 0 ) THEN
  PNSDERS => ZGTF(2*KST_nsders-1:2*(KST_nsders+KF_SCDERS-1),:)
  PEWDERS => ZGTF(2*KST_ewders-1:2*(KST_ewders+KF_SCDERS-1),:)
ENDIF
IF (LUVDER) THEN
  PUVDERS => ZGTF(2*KST_uvders-1:2*(KST_uvders+2*KF_UV-1),:)
ENDIF

!$ACC DATA PRESENT(ZGTF) &
!$ACC& PRESENT(PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

DO KGL=IBEG,IEND,IINC
  
IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = F%RACTHE(IGLG)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2  = F%RACTHE(IGLG)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRBT*ASIN(1.0_JPRBT)
  ZACHTE2 = 1._JPRBT
  ZACHTE  = F%RACTHE2(IGLG)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRBT)

    !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IR,II,ZAMP,ZPHASE)
    DO JF=1,KF_SCALARS
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        
        ! calculate amplitude and add phase shift then reconstruct A,B
        ZAMP = SQRT(PSCALAR(JF,IR)**2 + PSCALAR(JF,II)**2)
        ZPHASE = ATAN2(PSCALAR(JF,II),PSCALAR(JF,IR)) + REAL(JM,JPRBT)*ZSHIFT
        
        PSCALAR(2*JF-1,IR) =  ZAMP*COS(ZPHASE)
        PSCALAR(2*JF,  Ir) = ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
    IF(KF_SCDERS > 0)THEN
      !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IR,II,ZAMP,ZPHASE)
      DO JF=1,KF_SCALARS
        DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(PNSDERS(JF,IR)**2 + PNSDERS(JF,II)**2)
          ZPHASE = ATAN2(PNSDERS(JF,II),PNSDERS(JF,IR)) + REAL(JM,JPRBT)*ZSHIFT
          PNSDERS(2*JF-1,IR) =  ZAMP*COS(ZPHASE)
          PNSDERS(2*JF,  Ir) = ZAMP*SIN(ZPHASE)
        ENDDO
      ENDDO
    ENDIF
    !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IR,II,ZAMP,ZPHASE)
    DO JF=1,2*KF_UV
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        ! calculate amplitude and phase shift and reconstruct A,B
        ZAMP = SQRT(PUV(JF,IR)**2 + PUV(JF,II)**2)
        ZPHASE = ATAN2(PUV(JF,II),PUV(JF,IR)) + REAL(JM,JPRBT)*ZSHIFT
        PUV(2*JF-1,IR) =  ZAMP*COS(ZPHASE)
        PUV(2*JF,  Ir) =  ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
  !$ACC PARALLEL LOOP DEFAULT(NONE)
  DO JLON=ISTAGTF+1,ISTAGTF+2*IMEN+1
    DO JF=1,2*KF_UV
      PUV(2*JF-1,JLON) = PUV(2*JF-1,JLON)*ZACHTE2
      PUV(2*JF,  JLON) = PUV(2*JF  ,JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
  !$ACC PARALLEL LOOP DEFAULT(NONE)
  DO JLON=ISTAGTF+1,ISTAGTF+2*IMEN+1
    DO JF=1,KF_SCALARS
      PNSDERS(2*JF-1,JLON) = PNSDERS(2*JF-1,JLON)*ZACHTE2
      PNSDERS(2*JF,  JLON) = PNSDERS(2*JF,  JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  !$ACC PARALLEL LOOP PRIVATE(IR) DEFAULT(NONE)
  DO JM=0,IMEN
    DO JF=1,2*KF_UV
      IR = ISTAGTF+2*JM+1
      PUVDERS(2*JF-1,IR) = -PUV(2*JF,IR)*ZACHTE2*REAL(JM,JPRBT)
      PUVDERS(2*JF,  IR) =  PUV(2*JF-1,IR)*ZACHTE2*REAL(JM,JPRBT)
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
  !$ACC PARALLEL LOOP PRIVATE(IR) DEFAULT(NONE)
  DO JM=0,IMEN
    DO JF=1,KF_SCALARS
      IR = ISTAGTF+2*JM+1
      PEWDERS(2*JF-1,IR) = -PSCALAR(2*JF,IR)*ZACHTE2*REAL(JM,JPRBT)
      PEWDERS(2*JF,  IR) =  PSCALAR(2*JF-1,IR)*ZACHTE2*REAL(JM,JPRBT)
    ENDDO
  ENDDO
ENDIF

enddo
!$ACC END DATA
!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
