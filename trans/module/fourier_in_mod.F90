MODULE FOURIER_IN_MOD
CONTAINS
SUBROUTINE FOURIER_IN(PREEL,KFIELDS)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

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

#include "tsmbkind.h"

USE TPM_DISTR
USE TPM_TRANS
USE TPM_GEOMETRY

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KFIELDS

REAL_B, INTENT(OUT) :: PREEL(:,:)

INTEGER_M :: JGL,JM,JF,IGLG,IPROC,IR,II,ISTA

!     ------------------------------------------------------------------

!$OMP PARALLEL PRIVATE(JGL,IGLG,JM,IPROC,IR,II,ISTA,JF)
!$OMP DO SCHEDULE(STATIC,1)
DO JGL=1,D%NDGL_FS
  IGLG = D%NPTRLS(MYSETW)+JGL-1
  DO JM=0,G%NMEN(IGLG)
    IPROC = D%NPROCM(JM)
    IR    = 2*JM+1+D%NSTAGTF(JGL)
    II    = 2*JM+2+D%NSTAGTF(JGL)
    ISTA  = (D%NSTAGT0B(D%MSTABF(IPROC))+D%NPNTGTB0(JM,JGL))*2*KFIELDS
    DO JF=1,KFIELDS
      PREEL(JF,IR) = FOUBUF(ISTA+2*JF-1)
      PREEL(JF,II) = FOUBUF(ISTA+2*JF  )
    ENDDO
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_IN
END MODULE FOURIER_IN_MOD

