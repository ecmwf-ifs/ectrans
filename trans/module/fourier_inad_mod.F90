MODULE FOURIER_INAD_MOD
CONTAINS
SUBROUTINE FOURIER_INAD(PREEL,KFIELDS)

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

#include "tsmbkind.h"

USE TPM_DISTR
USE TPM_TRANS
USE TPM_GEOMETRY

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KFIELDS

REAL_B, INTENT(IN) :: PREEL(:,:)

INTEGER_M :: JGL,JM,JF,IGLG,IPROC,IR,II,ISTA

!     ------------------------------------------------------------------

#ifndef HLOMP
!$OMP PARALLEL PRIVATE(JGL,IGLG,JM,IPROC,IR,II,ISTA,JF)
!$OMP DO SCHEDULE(STATIC,1)
#endif
DO JGL=1,D%NDGL_FS
  IGLG = D%NPTRLS(MYSETW)+JGL-1
  DO JM=0,G%NMEN(IGLG)
    IPROC = D%NPROCM(JM)
    IR    = 2*JM+1+D%NSTAGTF(JGL)
    II    = 2*JM+2+D%NSTAGTF(JGL)
    ISTA  = (D%NSTAGT0B(D%MSTABF(IPROC))+D%NPNTGTB0(JM,JGL))*2*KFIELDS
    DO JF=1,KFIELDS
      FOUBUF(ISTA+2*JF-1) = PREEL(JF,IR)
      FOUBUF(ISTA+2*JF  ) = PREEL(JF,II)
    ENDDO
  ENDDO
ENDDO
#ifndef HLOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_INAD
END MODULE FOURIER_INAD_MOD

