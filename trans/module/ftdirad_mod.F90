MODULE FTDIRAD_MOD
CONTAINS
SUBROUTINE FTDIRAD(PREEL,KFIELDS,KGL)


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

INTEGER_M,INTENT(IN)  :: KFIELDS,KGL
REAL_B, INTENT(INOUT) :: PREEL(:,:)

INTEGER_M :: IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN,IOFF
REAL_B :: ZMUL
!     ------------------------------------------------------------------

IJUMP = 1

IGLG = D%NPTRLS(MYSETW)+KGL-1
ILOEN = G%NLOEN(IGLG)
IST  = 2*(G%NMEN(IGLG)+1)+1
ILEN = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+IOFF+JJ-1) = _ZERO_
  ENDDO
ENDDO

CALL FFT992(PREEL(1,IOFF+1),T%TRIGS(1,KGL),&
 &T%NFAX(1,KGL),KFIELDS,IJUMP,ILOEN,KFIELDS,1)

  ! Change of metric (not in forward routine)

ZMUL = _ONE_/ILOEN
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF+JJ) = PREEL(JF,IOFF+JJ)*ZMUL
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIRAD
END MODULE FTDIRAD_MOD
