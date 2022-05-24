! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
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

USE PARKIND_ECTRANS ,ONLY : JPIM, JPRBT
USE TPM_DISTR       ,ONLY : D,D_NSTAGTF,D_NPTRLS, MYSETW,  MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G, G_NLOEN, G_NMEN
USE TPM_GEN         ,ONLY : NOUT, LSYNC_TRANS
USE TPM_HICFFT      ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT, EXECUTE_INV_FFT
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE TPM_STATS       ,ONLY : GSTATS => GSTATS_NVTX
USE DEVICE_MOD
USE ISO_C_BINDING

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_REAL(:)
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_COMPLEX(:)

INTEGER(KIND=JPIM) :: IGLG,KGL,ISTAT

!     ------------------------------------------------------------------

ALLOCATE(PREEL_REAL(KFIELD*D%NLENGTF))
#ifdef ACCGPU
!$ACC ENTER DATA CREATE(PREEL_REAL)
!$ACC DATA PRESENT(PREEL_REAL,PREEL_COMPLEX)
#endif

IF (LSYNC_TRANS) THEN
  CALL GSTATS(440,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(440,1)
ENDIF
CALL GSTATS(423,0)
CALL EXECUTE_INV_FFT(PREEL_COMPLEX(:),PREEL_REAL(:),KFIELD, &
    & LOENS=G%NLOEN(D%NPTRLS(MYSETW):D%NPTRLS(MYSETW)+D%NDGL_FS-1), &
    & OFFSETS=D%NSTAGTF(1:D%NDGL_FS))

ISTAT = DEVICE_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL GSTATS(443,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(443,1)
ENDIF
CALL GSTATS(423,1)

!     ------------------------------------------------------------------
!

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif

#ifdef OMPGPU
#endif
#ifdef ACCGPU
!$ACC EXIT DATA DELETE(PREEL_COMPLEX)
#endif
DEALLOCATE(PREEL_COMPLEX)


END SUBROUTINE FTINV
END MODULE FTINV_MOD
