! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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

!     Externals.  FFTW - FFT routine
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
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW 
!        R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_GEOMETRY    ,ONLY : G
#ifdef WITH_FFT992
USE TPM_FFT992      ,ONLY : T992
USE BLUESTEIN_MOD   ,ONLY : BLUESTEIN_FFT
#endif
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_DIM         ,ONLY : R
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,JJ,JF,ILOEN
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN,ITYPE
REAL(KIND=JPRB) :: ZMUL
!     ------------------------------------------------------------------

ITYPE = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
IST   = 2*(G%NMEN(IGLG)+1)+1
ILOEN = G%NLOEN(IGLG)+R%NNOEXTZL
ILEN  = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)+1
IRLEN = ILOEN
ICLEN = (IRLEN/2+1)*2

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+IOFF-1+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

#ifdef WITH_FFTW
CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,TW%LALL_FFTW,PREEL)
#elif defined(WITH_FFT992)
IF( T992%LFFT992 .AND. T992%LUSEFFT992(KGL) ) THEN
  CALL FFT992(PREEL(1,IOFF),T992%TRIGS(1,KGL),T992%NFAX(1,KGL),KFIELDS,1,IRLEN,KFIELDS,ITYPE)
ELSEIF( T992%LFFT992 ) THEN
  CALL BLUESTEIN_FFT(T992%TB,IRLEN,ITYPE,KFIELDS,PREEL(1:KFIELDS,IOFF:IOFF+ICLEN-1))
ELSE
  CALL ABORT_TRANS('FTDIRAD: FFT992 fallback requires FFTW support in this build')
ENDIF
#else
CALL ABORT_TRANS('FTDIRAD: no Fourier backend available in this build')
#endif

  ! Change of metric (not in forward routine)

ZMUL = 1.0_JPRB/ILOEN
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF-1+JJ) = PREEL(JF,IOFF-1+JJ)*ZMUL
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIRAD
END MODULE FTDIRAD_MOD
