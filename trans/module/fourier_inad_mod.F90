MODULE FOURIER_INAD_MOD
CONTAINS
SUBROUTINE FOURIER_INAD(PREEL,KFIELDS,KGL)

!**** *FOURIER_INAD* - Copy fourier data from buffer to local array - adjoint

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_INAD(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR
USE TPM_TRANS
USE TPM_GEOMETRY

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL

REAL(KIND=JPRB), INTENT(IN) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,IR,II,ISTA

!     ------------------------------------------------------------------

IGLG = D%NPTRLS(MYSETW)+KGL-1
DO JM=0,G%NMEN(IGLG)
  IPROC = D%NPROCM(JM)
  IR    = 2*JM+1+D%NSTAGTF(KGL)
  II    = 2*JM+2+D%NSTAGTF(KGL)
  ISTA  = (D%NSTAGT0B(D%MSTABF(IPROC))+D%NPNTGTB0(JM,KGL))*2*KFIELDS
  DO JF=1,KFIELDS
    FOUBUF(ISTA+2*JF-1) = PREEL(JF,IR)
    FOUBUF(ISTA+2*JF  ) = PREEL(JF,II)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_INAD
END MODULE FOURIER_INAD_MOD

