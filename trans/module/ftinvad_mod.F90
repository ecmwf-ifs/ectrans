MODULE FTINVAD_MOD
CONTAINS
SUBROUTINE FTINVAD(PREEL,KFIELDS,KGL)


!**** *FTINVAD - Inverse Fourier transform - adjoint

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINVAD(..)

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
!        R. El Khatib 01-Sep-2015 Better modularity for FFTW
!     ------------------------------------------------------------------

USE, INTRINSIC :: ISO_C_BINDING

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_DIM         ,ONLY : R
!USE TPMALD_DIM      ,ONLY : RALD
!USE TPM_TRANS
USE TPM_GEOMETRY    ,ONLY : G
!USE EXTPER_MOD      ,ONLY : EXTPER
USE TPM_FFT         ,ONLY : T
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_DIM         ,ONLY : R
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

#include "fftw3.f03.h"

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(OUT)  :: PREEL(:,:)

REAL(KIND=JPRB),ALLOCATABLE :: ZFFT(:,:)
INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
INTEGER(KIND=JPIB) :: IPLAN_R2C

!     ------------------------------------------------------------------

ITYPE=-1
IJUMP = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
ILOEN = G%NLOEN(IGLG)+R%NNOEXTZL
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)+1

  ! Change of metric (not in forward routine)
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF+JJ-1) = PREEL(JF,IOFF+JJ-1)*ILOEN
  ENDDO
ENDDO
  
IF( T%LFFT992 )THEN
  CALL FFT992(PREEL(1,IOFF),T%TRIGS(1,KGL),&
   &T%NFAX(1,KGL),KFIELDS,IJUMP,ILOEN,KFIELDS,ITYPE)
#ifdef WITH_FFTW
ELSEIF( TW%LFFTW )THEN
    IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
    ICLEN=(IRLEN/2+1)*2
    CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,PREEL)
#endif
ELSE
  CALL ABORT_TRANS('FTINVAD: NO FFT PACKAGE SELECTED')
ENDIF

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST-1+IOFF+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTINVAD
END MODULE FTINVAD_MOD
