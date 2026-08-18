! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_TRANSFORM_MOD

USE PARKIND1, ONLY: JPIM

PRIVATE
PUBLIC :: FOURIER_TRANSFORM, JP_DIRECT, JP_INVERSE

INTEGER(KIND=JPIM), PARAMETER :: JP_DIRECT  = -1
INTEGER(KIND=JPIM), PARAMETER :: JP_INVERSE =  1

CONTAINS

SUBROUTINE FOURIER_TRANSFORM(KTYPE, PREEL, KFIELDS, KGL)

!**** *FOURIER_TRANSFORM - Fourier transform

!     Purpose. Routine for Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FOURIER_TRANSFORM(..)

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
!        G. Radnoti 01-04-24 : 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!        R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     ------------------------------------------------------------------

USE PARKIND1,     ONLY: JPIM, JPRB
USE TPM_DISTR,    ONLY: D, MYSETW
USE TPM_GEOMETRY, ONLY: G
USE TPM_FFTW,     ONLY: TW, EXEC_FFTW
USE TPM_DIM,      ONLY: R

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KTYPE
INTEGER(KIND=JPIM), INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM), INTENT(IN)    :: KGL
REAL(KIND=JPRB),    INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG, IST, ILEN, JJ, JF, IST1, IOFF, IRLEN, ICLEN

!     ------------------------------------------------------------------

! Note: R%NNOEXTZL only has relevance for limited-area transforms
! It is zero for global transforms

! Complex coefficients beyond the NMEN cutoff are zeroed before or after the transform, depending
! on the direction.

IGLG  = D%NPTRLS(MYSETW) + KGL - 1 ! Global latitude index
IST   = 2 * (G%NMEN(IGLG) + 1) + 1 ! Starting index of first complex coefficient to be zeroed
ILEN  = G%NLOEN(IGLG) + R%NNOEXTZL + 3 - IST ! Number of complex coefficients to be zeroed
IST1 = 1 ! Starting index of first complex coefficient to be zeroed, without offset
IF (G%NLOEN(IGLG) == 1) IST1 = 0

! If inverse, zero the complex coefficients beyond the NMEN cutoff before transforming
IF (KTYPE == 1) THEN
  DO JJ = IST1, ILEN
    PREEL(1:KFIELDS, IST + D%NSTAGTF(KGL) + JJ - 1) = 0.0_JPRB
  ENDDO
ENDIF

IF (G%NLOEN(IGLG) > 1) THEN
  IOFF = D%NSTAGTF(KGL) + 1 ! Offset to this latitude in packed PREEL array
  IRLEN = G%NLOEN(IGLG) + R%NNOEXTZL ! Real length of this FFT
  ICLEN = (IRLEN / 2 + 1) * 2 ! Complex length of this FFT

  CALL EXEC_FFTW(KTYPE, IRLEN, ICLEN, IOFF, KFIELDS, TW%LALL_FFTW, PREEL)
ENDIF

! If direct, zero the complex coefficients beyond the NMEN cutoff after transforming
IF (KTYPE == -1) THEN
  DO JJ = IST1, ILEN
    PREEL(1:KFIELDS, IST + D%NSTAGTF(KGL) + JJ - 1) = 0.0_JPRB
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_TRANSFORM
END MODULE FOURIER_TRANSFORM_MOD
