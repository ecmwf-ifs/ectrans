MODULE FOURIER_OUT_MOD
CONTAINS
SUBROUTINE FOURIER_OUT(PREEL,KFIELDS)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

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

REAL_B,   INTENT(IN) :: PREEL(:,:)
INTEGER_M,INTENT(IN) :: KFIELDS

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
    ISTA  = (D%NSTAGT1B(D%MSTABF(IPROC))+D%NPNTGTB0(JM,JGL))*2*KFIELDS
    DO JF=1,KFIELDS
      FOUBUF_IN(ISTA+2*JF-1) = PREEL(JF,IR)
      FOUBUF_IN(ISTA+2*JF  ) = PREEL(JF,II)
    ENDDO
  ENDDO
ENDDO
#ifndef HLOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_OUT
END MODULE FOURIER_OUT_MOD

