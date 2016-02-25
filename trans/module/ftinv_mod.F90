MODULE FTINV_MOD
CONTAINS
SUBROUTINE FTINV(PREEL,KFIELDS,KGL)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

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
!        G. Radnoti 01-04-24 : 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
!USE TPM_TRANS
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFT         ,ONLY : T
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_DIM         ,ONLY : R
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(OUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

!     ------------------------------------------------------------------

ITYPE = 1
IJUMP = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = G%NLOEN(IGLG)+R%NNOEXTZL+3-IST
IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0

DO JJ=IST1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

IF (G%NLOEN(IGLG)>1) THEN
  IOFF=D%NSTAGTF(KGL)+1
  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  IF( T%LFFT992 )THEN
    CALL FFT992(PREEL(1,IOFF),T%TRIGS(1,KGL),&
     &T%NFAX(1,KGL),KFIELDS,IJUMP,IRLEN,KFIELDS,ITYPE)
#ifdef WITH_FFTW
  ELSEIF( TW%LFFTW )THEN
    ICLEN=(IRLEN/2+1)*2
    CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,LL_ALL,PREEL)
#endif
  ELSE
    CALL ABORT_TRANS('FTINV: NO FFT PACKAGE SELECTED')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTINV
END MODULE FTINV_MOD
