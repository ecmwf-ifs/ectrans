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
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
USE TPM_DIM         ,ONLY : R
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS,KGL
REAL(KIND=JPRB), POINTER, CONTIGUOUS, INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,JJ,IST1
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

  CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,TW%LALL_FFTW,PREEL)
ENDIF

IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0
DO JJ=IST1,ILEN
  PREEL(1:KFIELDS,IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
