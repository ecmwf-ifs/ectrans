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
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms

!     ------------------------------------------------------------------

USE, INTRINSIC :: ISO_C_BINDING

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
!USE TPM_TRANS
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFT         ,ONLY : T
USE TPM_FFTW        ,ONLY : TW, CREATE_PLAN_FFTW
USE TPM_DIM         ,ONLY : R
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE

#include "fftw3.f03.h"

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)

REAL(KIND=JPRB),ALLOCATABLE :: ZFFT(:,:)
INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN
INTEGER(KIND=JPIB) :: IPLAN_C2R
REAL(KIND=JPRB) :: ZMUL
!     ------------------------------------------------------------------

IJUMP = 1
IGLG = D%NPTRLS(MYSETW)+KGL-1
ILOEN = G%NLOEN(IGLG)+R%NNOEXTZL
IST  = 2*(G%NMEN(IGLG)+1)+1
ILEN = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+IOFF+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

IF( T%LFFT992 )THEN
  CALL FFT992(PREEL(1,IOFF+1),T%TRIGS(1,KGL),&
   &T%NFAX(1,KGL),KFIELDS,IJUMP,ILOEN,KFIELDS,1)
ELSEIF( TW%LFFTW )THEN
  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  ICLEN=(IRLEN/2+1)*2
  CALL CREATE_PLAN_FFTW(IPLAN_C2R,1,IRLEN,KFIELDS)
  ALLOCATE(ZFFT(IOFF+1:IOFF+ICLEN,KFIELDS))
  DO JF=1,KFIELDS
    DO JJ=1,ICLEN
      ZFFT(IOFF+JJ,JF) =PREEL(JF,IOFF+JJ)
    ENDDO
  ENDDO
  CALL DFFTW_EXECUTE_DFT_C2R(IPLAN_C2R,ZFFT,ZFFT)
  DO JJ=1,IRLEN
    DO JF=1,KFIELDS
      PREEL(JF,IOFF+JJ)=ZFFT(IOFF+JJ,JF)
    ENDDO
  ENDDO
  DEALLOCATE(ZFFT)
ELSE
  CALL ABORT_TRANS('FTDIRAD: NO FFT PACKAGE SELECTED')
ENDIF

  ! Change of metric (not in forward routine)

ZMUL = 1.0_JPRB/ILOEN
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF+JJ) = PREEL(JF,IOFF+JJ)*ZMUL
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIRAD
END MODULE FTDIRAD_MOD
