! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
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

USE TPM_TRANS       ,ONLY : LUVDER, LATLON, ZGTF, LVORGP, LDIVGP
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

REAL(KIND=JPRBT) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRBT) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JF,IGLG,II,IR,JM

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
!DEBUGGING:
integer :: i,J,maxi,maxj,ist
real :: maxv

INTEGER(JPIM) :: ZGTF_START(8)
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_VOR = 1
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_DIV = 2
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_UV = 3
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_SCALAR = 4
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_NSDERS = 5
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_UVDERS = 6
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_EWDERS = 7
INTEGER(JPIM), PARAMETER :: ZGTF_START_INDEX_END = 8
INTEGER(JPIM) :: JF_UV, JF_SCALAR

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

IST = 1
ZGTF_START(ZGTF_START_INDEX_VOR) = IST
IF (LVORGP) THEN
  IST = IST+KF_UV
ENDIF
ZGTF_START(ZGTF_START_INDEX_DIV) = IST
IF (LDIVGP) THEN
  IST = IST+KF_UV
ENDIF
ZGTF_START(ZGTF_START_INDEX_UV) = IST
IST = IST+2*KF_UV
ZGTF_START(ZGTF_START_INDEX_SCALAR) = IST
IST = IST+KF_SCALARS
ZGTF_START(ZGTF_START_INDEX_NSDERS) = IST
IST = IST+KF_SCDERS
ZGTF_START(ZGTF_START_INDEX_UVDERS) = IST
IF (LUVDER) THEN
  IST = IST+2*KF_UV
ENDIF
ZGTF_START(ZGTF_START_INDEX_EWDERS) = IST
IF (KF_SCDERS > 0) THEN
  IST = IST+KF_SCDERS
ENDIF
ZGTF_START(ZGTF_START_INDEX_END) = IST

!$ACC DATA PRESENT(ZGTF,D,G,F) COPYIN(ZGTF_START)

IF( LATLON.AND.S%LDLL.AND.S%LSHIFTLL ) THEN
  PRINT *, "This is not implemented yet! LATLON.AND.S%LDLL.AND.S%LSHIFTLL"
  STOP 128 ! not implemented
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

!*       1.1      U AND V.

IF (ZGTF_START(ZGTF_START_INDEX_UV) /= ZGTF_START(ZGTF_START_INDEX_UV+1)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=ZGTF_START(ZGTF_START_INDEX_UV),ZGTF_START(ZGTF_START_INDEX_UV+1)-1

      IGLG    = D%NPTRLS(MYSETW)+KGL-1
      IMEN    = G%NMEN(IGLG)
      ISTAGTF = D%NSTAGTF(KGL)
      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,2*IMEN
        IR = ISTAGTF+JM+1
        ZGTF(2*JF-1,IR) = ZGTF(2*JF-1,IR)*ZACHTE2
        ZGTF(2*JF,  IR) = ZGTF(2*JF  ,IR)*ZACHTE2
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF (ZGTF_START(ZGTF_START_INDEX_NSDERS) /= ZGTF_START(ZGTF_START_INDEX_NSDERS+1)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=ZGTF_START(ZGTF_START_INDEX_NSDERS),ZGTF_START(ZGTF_START_INDEX_NSDERS+1)-1

      IGLG    = D%NPTRLS(MYSETW)+KGL-1
      IMEN    = G%NMEN(IGLG)
      ISTAGTF = D%NSTAGTF(KGL)
      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,2*IMEN
        IR = ISTAGTF+JM+1
        ZGTF(2*JF-1,IR) = ZGTF(2*JF-1,IR)*ZACHTE2
        ZGTF(2*JF,  IR) = ZGTF(2*JF,  IR)*ZACHTE2
      ENDDO
    ENDDO
  ENDDO
ENDIF

    !     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF (ZGTF_START(ZGTF_START_INDEX_UVDERS) /= ZGTF_START(ZGTF_START_INDEX_UVDERS+1)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=ZGTF_START(ZGTF_START_INDEX_UVDERS),ZGTF_START(ZGTF_START_INDEX_UVDERS+1)-1

      IGLG    = D%NPTRLS(MYSETW)+KGL-1
      IMEN    = G%NMEN(IGLG)
      ISTAGTF = D%NSTAGTF(KGL)
      ZACHTE2  = F%RACTHE(IGLG)
      JF_UV = JF - ZGTF_START(ZGTF_START_INDEX_UVDERS) + ZGTF_START(ZGTF_START_INDEX_UV)

      !$ACC LOOP SEQ
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        ZGTF(2*JF-1,IR) = -ZGTF(2*JF_UV,IR)*ZACHTE2*REAL(JM,JPRBT)
        ZGTF(2*JF,  IR) =  ZGTF(2*JF_UV-1,IR)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES
IF (ZGTF_START(ZGTF_START_INDEX_EWDERS) /= ZGTF_START(ZGTF_START_INDEX_EWDERS+1)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=ZGTF_START(ZGTF_START_INDEX_EWDERS),ZGTF_START(ZGTF_START_INDEX_EWDERS+1)-1

      IGLG    = D%NPTRLS(MYSETW)+KGL-1
      IMEN    = G%NMEN(IGLG)
      ISTAGTF = D%NSTAGTF(KGL)
      ZACHTE2  = F%RACTHE(IGLG)
      JF_SCALAR = JF - ZGTF_START(ZGTF_START_INDEX_EWDERS) + ZGTF_START(ZGTF_START_INDEX_SCALAR)

      !$ACC LOOP SEQ
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        ZGTF(2*JF-1,IR) = -ZGTF(2*JF_SCALAR,IR)*ZACHTE2*REAL(JM,JPRBT)
        ZGTF(2*JF,  IR) =  ZGTF(2*JF_SCALAR-1,IR)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!$ACC WAIT(1)

!$ACC END DATA
!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
