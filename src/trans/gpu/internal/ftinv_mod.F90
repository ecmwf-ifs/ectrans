! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV_MOD
CONTAINS
SUBROUTINE FTINV(PREEL_COMPLEX,PREEL_REAL,KFIELD)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELD   - number of fields

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
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!     ------------------------------------------------------------------

USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE PARKIND_ECTRANS ,ONLY : JPIM, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFTC        ,ONLY : EXECUTE_INV_FFT
USE CUDA_DEVICE_MOD
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_REAL(:)
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_COMPLEX(:)

INTEGER(KIND=JPIM) :: IGLG,KGL,IRET

!     ------------------------------------------------------------------

PREEL_REAL => PREEL_COMPLEX

!$ACC DATA PRESENT(PREEL_REAL,PREEL_COMPLEX)

IF (LSYNC_TRANS) THEN
  CALL GSTATS(440,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(440,1)
ENDIF
CALL GSTATS(423,0)
CALL EXECUTE_INV_FFT(PREEL_COMPLEX(:),PREEL_REAL(:),KFIELD, &
    & LOENS=G%NLOEN(D%NPTRLS(MYSETW):D%NPTRLS(MYSETW)+D%NDGL_FS-1), &
    & OFFSETS=D%NSTAGTF(1:D%NDGL_FS))
IRET = CUDA_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL GSTATS(443,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(443,1)
ENDIF
CALL GSTATS(423,1)

!$ACC END DATA

NULLIFY(PREEL_COMPLEX)

!     ------------------------------------------------------------------

END SUBROUTINE FTINV
END MODULE FTINV_MOD
