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
SUBROUTINE FSC(KF_UV, KF_SCALARS, PUV, PSCALARS, PSCALARS_NSDER, PUV_EWDER, PSCALARS_EWDER)

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
USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC, D_NSTAGTF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN
USE TPM_FIELDS      ,ONLY : F
USE TPM_FLT                ,ONLY: S
USE TPM_GEN, ONLY: NOUT

USE TPM_TRANS       ,ONLY : LATLON
!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KGL

REAL(KIND=JPRBT) :: ZACHTE2
INTEGER(KIND=JPIM) :: IOFF,OFFSET_VAR
INTEGER(KIND=JPIM) :: JF,IGLG,JM
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV, KF_SCALARS
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: &
    & PUV(:,:), PSCALARS(:,:), PSCALARS_NSDER(:,:), PUV_EWDER(:,:), PSCALARS_EWDER(:,:)

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

!$ACC DATA PRESENT(D,G,F,D_NSTAGTF,G_NMEN,PUV,PSCALARS,PSCALARS_NSDER,PUV_EWDER,PSCALARS_EWDER)

IF( LATLON.AND.S%LDLL.AND.S%LSHIFTLL ) THEN
  PRINT *, "This is not implemented yet! LATLON.AND.S%LDLL.AND.S%LSHIFTLL"
  STOP 128 ! not implemented
ENDIF
 
OFFSET_VAR=D%NPTRLS(MYSETW)

!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

!*       1.1      U AND V.

!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
DO KGL=IBEG,IEND,IINC
  DO JF=1,2*KF_UV
    IGLG    = OFFSET_VAR+KGL-1
    IOFF = D_NSTAGTF(KGL)+1

    ZACHTE2  = F%RACTHE(IGLG)

    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      PUV(2*JF-1,2*JM+IOFF) = PUV(2*JF-1,2*JM+IOFF)*ZACHTE2
      PUV(2*JF,  2*JM+IOFF) = PUV(2*JF  ,2*JM+IOFF)*ZACHTE2
    ENDDO
  ENDDO
ENDDO

!*      1.2      N-S DERIVATIVES

IF (ASSOCIATED(PSCALARS_NSDER)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=1,KF_SCALARS
      IGLG = OFFSET_VAR+KGL-1
      IOFF = D_NSTAGTF(KGL)+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        PSCALARS_NSDER(2*JF-1,2*JM+IOFF) = PSCALARS_NSDER(2*JF-1,2*JM+IOFF)*ZACHTE2
        PSCALARS_NSDER(2*JF,  2*JM+IOFF) = PSCALARS_NSDER(2*JF,  2*JM+IOFF)*ZACHTE2
      ENDDO
    ENDDO
  ENDDO
ENDIF

    !     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF (ASSOCIATED(PUV_EWDER)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=1,2*KF_UV
      IGLG = OFFSET_VAR+KGL-1
      IOFF = D_NSTAGTF(KGL)+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        PUV_EWDER(2*JF-1,2*JM+IOFF) = -PUV(2*JF,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
        PUV_EWDER(2*JF,  2*JM+IOFF) =  PUV(2*JF-1,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES
IF (ASSOCIATED(PSCALARS_EWDER)) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=1,KF_SCALARS
      IGLG = OFFSET_VAR+KGL-1
      IOFF = D_NSTAGTF(KGL)+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        PSCALARS_EWDER(2*JF-1,2*JM+IOFF) = -PSCALARS(2*JF,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
        PSCALARS_EWDER(2*JF,  2*JM+IOFF) =  PSCALARS(2*JF-1,2*JM+IOFF)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!$ACC WAIT(1)

!$ACC END DATA
!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
