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
USE TPM_FFTW        ,ONLY : TW, CREATE_PLAN_FFTW
#endif
USE TPM_DIM         ,ONLY : R
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

#include "fftw3.f03.h"

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(OUT)  :: PREEL(:,:)

REAL(KIND=JPRB), POINTER :: ZFFT(:,:)
REAL(KIND=JPRB), POINTER :: ZFFT1(:)
TYPE(C_PTR) :: ZFFTP, ZFFT1P
INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN
INTEGER(KIND=JPIB) :: IPLAN_R2C
INTEGER(KIND=JPIB) :: IPLAN_R2C1
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

!     ------------------------------------------------------------------

IJUMP = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
ILOEN = G%NLOEN(IGLG)+R%NNOEXTZL
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)

  ! Change of metric (not in forward routine)
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF+JJ) = PREEL(JF,IOFF+JJ)*ILOEN
  ENDDO
ENDDO
  
IF( T%LFFT992 )THEN
  CALL FFT992(PREEL(1,IOFF+1),T%TRIGS(1,KGL),&
   &T%NFAX(1,KGL),KFIELDS,IJUMP,ILOEN,KFIELDS,-1)
#ifdef WITH_FFTW
ELSEIF( TW%LFFTW )THEN
  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  ICLEN=(IRLEN/2+1)*2
  IF( LL_ALL )THEN
    CALL CREATE_PLAN_FFTW(IPLAN_R2C,-1,IRLEN,KFIELDS)
    ZFFTP=FFTW_ALLOC_COMPLEX(INT(ICLEN/2*KFIELDS,C_SIZE_T))
    CALL C_F_POINTER(ZFFTP,ZFFT,[ICLEN,KFIELDS])
    DO JF=1,KFIELDS
      DO JJ=1,IRLEN
        ZFFT(JJ,JF) =PREEL(JF,IOFF+JJ)
      ENDDO
    ENDDO
    CALL DFFTW_EXECUTE_DFT_R2C(IPLAN_R2C,ZFFT,ZFFT)
    DO JJ=1,ICLEN
      DO JF=1,KFIELDS
        PREEL(JF,IOFF+JJ)=ZFFT(JJ,JF)/REAL(IRLEN,JPRB)
      ENDDO
    ENDDO
    CALL FFTW_FREE(ZFFTP)
  ELSE
    CALL CREATE_PLAN_FFTW(IPLAN_R2C1,-1,IRLEN,1)
    ZFFT1P=FFTW_ALLOC_COMPLEX(INT(ICLEN/2,C_SIZE_T))
    CALL C_F_POINTER(ZFFT1P,ZFFT1,[ICLEN])
    DO JF=1,KFIELDS
      DO JJ=1,IRLEN
        ZFFT1(JJ) =PREEL(JF,IOFF+JJ)
      ENDDO
      CALL DFFTW_EXECUTE_DFT_R2C(IPLAN_R2C1,ZFFT1,ZFFT1)
      DO JJ=1,ICLEN
        PREEL(JF,IOFF+JJ)=ZFFT1(JJ)/REAL(IRLEN,JPRB)
      ENDDO
    ENDDO
    CALL FFTW_FREE(ZFFT1P)
  ENDIF
#endif
ELSE
  CALL ABORT_TRANS('FTINVAD: NO FFT PACKAGE SELECTED')
ENDIF

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+IOFF+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTINVAD
END MODULE FTINVAD_MOD
