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
SUBROUTINE FSC(PREEL_COMPLEX, KF_UV, KF_SCALARS, KUV_OFFSET, &
        & KSCALARS_OFFSET, KSCALARS_NSDER_OFFSET, KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET)

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
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN, G_NLOEN
USE TPM_FIELDS      ,ONLY : F
USE TPM_FLT                ,ONLY: S
USE TPM_GEN, ONLY: NOUT

USE TPM_TRANS       ,ONLY : LATLON
!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KGL

REAL(KIND=JPRBT) :: ZACHTE2
INTEGER(KIND=JPIM) :: IOFF_COMPLEX,OFFSET_VAR
INTEGER(KIND=JPIM) :: JF,IGLG,JM
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

REAL(KIND=JPRBT), INTENT(INOUT) :: PREEL_COMPLEX(:,:)
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV, KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KUV_OFFSET, KSCALARS_OFFSET, KSCALARS_NSDER_OFFSET, KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET

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

!$ACC DATA PRESENT(D,G,F,D_NSTAGTF,G_NMEN,G_NLOEN,PREEL_COMPLEX)

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
    IOFF_COMPLEX = D_NSTAGTF(KGL)/2+1

    ZACHTE2  = F%RACTHE(IGLG)

    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      PREEL_COMPLEX(KUV_OFFSET+2*JF-1,JM+IOFF_COMPLEX) = &
          & PREEL_COMPLEX(KUV_OFFSET+2*JF-1,JM+IOFF_COMPLEX)*ZACHTE2
      PREEL_COMPLEX(KUV_OFFSET+2*JF,  JM+IOFF_COMPLEX) = &
          & PREEL_COMPLEX(KUV_OFFSET+2*JF  ,JM+IOFF_COMPLEX)*ZACHTE2
    ENDDO
  ENDDO
ENDDO

!*      1.2      N-S DERIVATIVES

IF (KSCALARS_NSDER_OFFSET >= 0) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=1,KF_SCALARS
      IGLG = OFFSET_VAR+KGL-1
      IOFF_COMPLEX = D_NSTAGTF(KGL)/2+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        PREEL_COMPLEX(KSCALARS_NSDER_OFFSET+2*JF-1,JM+IOFF_COMPLEX) = &
            & PREEL_COMPLEX(KSCALARS_NSDER_OFFSET+2*JF-1,JM+IOFF_COMPLEX)*ZACHTE2
        PREEL_COMPLEX(KSCALARS_NSDER_OFFSET+2*JF,  JM+IOFF_COMPLEX) = &
            & PREEL_COMPLEX(KSCALARS_NSDER_OFFSET+2*JF,  JM+IOFF_COMPLEX)*ZACHTE2
      ENDDO
    ENDDO
  ENDDO
ENDIF

    !     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF (KUV_EWDER_OFFSET >= 0) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=1,2*KF_UV
      IGLG = OFFSET_VAR+KGL-1
      IOFF_COMPLEX = D_NSTAGTF(KGL)/2+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        PREEL_COMPLEX(KUV_EWDER_OFFSET+2*JF-1,JM+IOFF_COMPLEX) = &
            & -PREEL_COMPLEX(KUV_OFFSET+2*JF,  JM+IOFF_COMPLEX)*ZACHTE2*REAL(JM,JPRBT)
        PREEL_COMPLEX(KUV_EWDER_OFFSET+2*JF,  JM+IOFF_COMPLEX) =  &
            &  PREEL_COMPLEX(KUV_OFFSET+2*JF-1,JM+IOFF_COMPLEX)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
      !$ACC LOOP SEQ
      DO JM=G_NMEN(IGLG)+1,(G_NLOEN(IGLG)+4)/2-1
        PREEL_COMPLEX(KUV_EWDER_OFFSET+2*JF-1,JM+IOFF_COMPLEX) = 0._JPRBT
        PREEL_COMPLEX(KUV_EWDER_OFFSET+2*JF,  JM+IOFF_COMPLEX) = 0._JPRBT
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES
IF (KSCALARS_EWDER_OFFSET > 0) THEN
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) ASYNC(1)
  DO KGL=IBEG,IEND,IINC
    DO JF=1,KF_SCALARS
      IGLG = OFFSET_VAR+KGL-1
      IOFF_COMPLEX = D_NSTAGTF(KGL)/2+1

      ZACHTE2  = F%RACTHE(IGLG)

      !$ACC LOOP SEQ
      DO JM=0,G_NMEN(IGLG)
        PREEL_COMPLEX(KSCALARS_EWDER_OFFSET+2*JF-1,JM+IOFF_COMPLEX) = &
            & -PREEL_COMPLEX(KSCALARS_OFFSET+2*JF  ,JM+IOFF_COMPLEX)*ZACHTE2*REAL(JM,JPRBT)
        PREEL_COMPLEX(KSCALARS_EWDER_OFFSET+2*JF,  JM+IOFF_COMPLEX) = &
            &  PREEL_COMPLEX(KSCALARS_OFFSET+2*JF-1,JM+IOFF_COMPLEX)*ZACHTE2*REAL(JM,JPRBT)
      ENDDO
      !$ACC LOOP SEQ
      DO JM=G_NMEN(IGLG)+1,(G_NLOEN(IGLG)+4)/2-1
        PREEL_COMPLEX(KSCALARS_EWDER_OFFSET+2*JF-1,JM+IOFF_COMPLEX) = 0._JPRBT
        PREEL_COMPLEX(KSCALARS_EWDER_OFFSET+2*JF,  JM+IOFF_COMPLEX) = 0._JPRBT
      ENDDO
    ENDDO
  ENDDO
ENDIF

!$ACC WAIT(1)

!$ACC END DATA
!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
