MODULE FTDIRAD_MOD
CONTAINS
SUBROUTINE FTDIRAD(PREEL,KFIELDS)


!**** *FTDIRAD - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform - adjoint
!     --------

!**   Interface.
!     ----------
!        CALL FTDIRAD(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------  
!                 

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DISTR
USE TPM_TRANS
USE TPM_GEOMETRY
USE TPM_FFT

IMPLICIT NONE

INTEGER_M,INTENT(IN)  :: KFIELDS
REAL_B, INTENT(INOUT) :: PREEL(:,:)

INTEGER_M :: JGL,IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN,IOFF

!     ------------------------------------------------------------------

IJUMP = 1

!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JGL,IGLG,ILOEN,IST,ILEN,IOFF,JJ,JF)
DO JGL=1,D%NDGL_FS
  IGLG = D%NPTRLS(MYSETW)+JGL-1
  ILOEN = G%NLOEN(IGLG)
  IST  = 2*(G%NMEN(IGLG)+1)+1
  ILEN = ILOEN+3-IST
  IOFF  = D%NSTAGTF(JGL)
  
  DO JJ=1,ILEN
    DO JF=1,KFIELDS
      PREEL(JF,IST+IOFF+JJ-1) = _ZERO_
    ENDDO
  ENDDO

  CALL FFT992(PREEL(1,IOFF+1),T%TRIGS(1,JGL),&
   &T%NFAX(1,JGL),KFIELDS,IJUMP,ILOEN,KFIELDS,1)

  ! Change of metric (not in forward routine)

  DO JJ=1,ILOEN
    DO JF=1,KFIELDS
      PREEL(JF,IOFF+JJ) = PREEL(JF,IOFF+JJ)/ILOEN
    ENDDO
  ENDDO
 

ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIRAD
END MODULE FTDIRAD_MOD
