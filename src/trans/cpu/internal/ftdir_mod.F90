! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(PREEL,KFIELDS,KGL)


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

!     Externals.  FFTW - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!        R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB

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
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS,KGL
REAL(KIND=JPRB), POINTER, CONTIGUOUS, INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE

!     ------------------------------------------------------------------

ITYPE=-1
IGLG = D%NPTRLS(MYSETW)+KGL-1
IST  = 2*(G%NMEN(IGLG)+1)+1
ILEN = G%NLOEN(IGLG)+R%NNOEXTZL+3-IST

IF (G%NLOEN(IGLG)>1) THEN 
  IOFF=D%NSTAGTF(KGL)+1
  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  ICLEN=(IRLEN/2+1)*2

#ifdef WITH_FFTW
  CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,TW%LALL_FFTW,PREEL)
#elif defined(WITH_FFT992)
  IF( T992%LFFT992 .AND. T992%LUSEFFT992(KGL) ) THEN
    CALL FFT992(PREEL(1,IOFF),T992%TRIGS(1,KGL),T992%NFAX(1,KGL),KFIELDS,1,IRLEN,KFIELDS,ITYPE)
  ELSEIF( T992%LFFT992 ) THEN
    CALL BLUESTEIN_FFT(T992%TB,IRLEN,ITYPE,KFIELDS,PREEL(1:KFIELDS,IOFF:IOFF+ICLEN-1))
    DO JJ=1,ICLEN
      DO JF=1,KFIELDS
        PREEL(JF,IOFF+JJ-1)=PREEL(JF,IOFF+JJ-1)/REAL(IRLEN,JPRB)
      ENDDO
    ENDDO
  ELSE
    CALL ABORT_TRANS('FTDIR: FFT992 fallback requires FFTW support in this build')
  ENDIF
#else
  CALL ABORT_TRANS('FTDIR: no Fourier backend available in this build')
#endif
ENDIF

IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0
DO JJ=IST1,ILEN
  PREEL(1:KFIELDS,IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
