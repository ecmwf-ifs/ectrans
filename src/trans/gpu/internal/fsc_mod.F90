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
SUBROUTINE FSC

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
USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC, D_NSTAGTF
USE TPM_FIELDS      ,ONLY : F, ZGTF_START, ZGTF_START_INDEX_UV, ZGTF_START_INDEX_NSDERS, &
    & ZGTF_START_INDEX_UVDERS, ZGTF_START_INDEX_EWDERS, ZGTF_START_INDEX_SCALAR
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN
USE TPM_FLT                ,ONLY: S
use tpm_gen, only: nout
!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KGL

REAL(KIND=JPRBT) :: ZACHTE2
INTEGER(KIND=JPIM) :: IOFF,OFFSET_VAR
INTEGER(KIND=JPIM) :: JF,IGLG,JM,JF_UV,JF_SCALAR
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

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

!$ACC DATA PRESENT(ZGTF,D,G,F,D_NSTAGTF,G_NMEN) COPYIN(ZGTF_START)

IF( LATLON.AND.S%LDLL.AND.S%LSHIFTLL ) THEN
  PRINT *, "This is not implemented yet! LATLON.AND.S%LDLL.AND.S%LSHIFTLL"
  STOP 128 ! not implemented
ENDIF
 
OFFSET_VAR=D%NPTRLS(MYSETW)

!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

!*       1.1      U AND V.

IF (ZGTF_START(ZGTF_START_INDEX_UV) /= ZGTF_START(ZGTF_START_INDEX_UV+1)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=ZGTF_START(ZGTF_START_INDEX_UV),ZGTF_START(ZGTF_START_INDEX_UV+1)-1
      IGLG    = OFFSET_VAR+KGL-1
      IOFF = D_NSTAGTF(KGL)+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        ZGTF(2*JF-1,2*JM+IOFF) = ZGTF(2*JF-1,2*JM+IOFF)*ZACHTE2
        ZGTF(2*JF,  2*JM+IOFF) = ZGTF(2*JF  ,2*JM+IOFF)*ZACHTE2
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF (ZGTF_START(ZGTF_START_INDEX_NSDERS) /= ZGTF_START(ZGTF_START_INDEX_NSDERS+1)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=ZGTF_START(ZGTF_START_INDEX_NSDERS),ZGTF_START(ZGTF_START_INDEX_NSDERS+1)-1
      IGLG = OFFSET_VAR+KGL-1
      IOFF = D_NSTAGTF(KGL)+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        ZGTF(2*JF-1,2*JM+IOFF) = ZGTF(2*JF-1,2*JM+IOFF)*ZACHTE2
        ZGTF(2*JF,  2*JM+IOFF) = ZGTF(2*JF,  2*JM+IOFF)*ZACHTE2
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
      IGLG = OFFSET_VAR+KGL-1
      IOFF = D_NSTAGTF(KGL)+1

      ZACHTE2  = F%RACTHE(IGLG)
      JF_UV = JF - ZGTF_START(ZGTF_START_INDEX_UVDERS) + ZGTF_START(ZGTF_START_INDEX_UV)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        ZGTF(2*JF-1,2*JM+IOFF) = -ZGTF(2*JF_UV,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
        ZGTF(2*JF,  2*JM+IOFF) =  ZGTF(2*JF_UV-1,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES
IF (ZGTF_START(ZGTF_START_INDEX_EWDERS) /= ZGTF_START(ZGTF_START_INDEX_EWDERS+1)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=ZGTF_START(ZGTF_START_INDEX_EWDERS),ZGTF_START(ZGTF_START_INDEX_EWDERS+1)-1
      IGLG = OFFSET_VAR+KGL-1
      IOFF = D_NSTAGTF(KGL)+1

      ZACHTE2  = F%RACTHE(IGLG)
      JF_SCALAR = JF - ZGTF_START(ZGTF_START_INDEX_EWDERS) + ZGTF_START(ZGTF_START_INDEX_SCALAR)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        ZGTF(2*JF-1,2*JM+IOFF) = -ZGTF(2*JF_SCALAR,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
        ZGTF(2*JF,  2*JM+IOFF) =  ZGTF(2*JF_SCALAR-1,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!$ACC WAIT(1)

!$ACC END DATA
!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
