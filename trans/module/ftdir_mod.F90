MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(PREEL,KFIELDS)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

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

INTEGER_M :: JGL,IGLG,IST,ILEN,IJUMP,JJ,JF

!     ------------------------------------------------------------------

IJUMP = 1

!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JGL,IGLG,IST,ILEN,JJ,JF)
DO JGL=1,D%NDGL_FS
  IGLG = D%NPTRLS(MYSETW)+JGL-1
  IST  = 2*(G%NMEN(IGLG)+1)+1
  ILEN = G%NLOEN(IGLG)+3-IST
  
  CALL FFT992(PREEL(1,D%NSTAGTF(JGL)+1),T%TRIGS(1,JGL),&
   &T%NFAX(1,JGL),NF_FS,IJUMP,G%NLOEN(IGLG),KFIELDS,-1)

  DO JJ=1,ILEN
    DO JF=1,KFIELDS
      PREEL(JF,IST+D%NSTAGTF(JGL)+JJ-1) = _ZERO_
    ENDDO
  ENDDO

ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
