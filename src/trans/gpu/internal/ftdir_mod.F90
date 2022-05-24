! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(PREEL_REAL,PREEL_COMPLEX,KFIELD)

!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

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
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT, EXECUTE_DIR_FFT
USE CUDA_DEVICE_MOD
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_REAL(:)
REAL(KIND=JPRBT), INTENT(OUT), POINTER :: PREEL_COMPLEX(:)

INTEGER(KIND=JPIM) :: IGLG,IOFF_REAL,IOFF_COMPLEX,KGL,IRET
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

!     ------------------------------------------------------------------

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF


ALLOCATE(PREEL_COMPLEX(KFIELD*D%NLENGTF))
!$ACC ENTER DATA CREATE(PREEL_COMPLEX)

!$ACC DATA PRESENT(PREEL_REAL,PREEL_COMPLEX)

IF (LSYNC_TRANS) THEN
  CALL GSTATS(430,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(430,1)
ENDIF
CALL GSTATS(413,0)
CALL EXECUTE_DIR_FFT(PREEL_REAL(:),PREEL_COMPLEX(:),KFIELD, &
    & LOENS=G%NLOEN(D%NPTRLS(MYSETW):D%NPTRLS(MYSETW)+D%NDGL_FS-1), &
    & OFFSETS=D%NSTAGTF(1:D%NDGL_FS))

DO KGL=IBEG,IEND,IINC

  ! NSTAGTF gives us space for NLOEN+3 elements
  ! In reality, at this point we need space for at most NLOEN+2 elements
  ! (in case NLOEN is even, otherwise NLOEN+1, due to the R2C definition)
  IOFF_REAL=D%NSTAGTF(KGL)+1
  IOFF_COMPLEX=D%NSTAGTF(KGL)/2+1
  IGLG = D%NPTRLS(MYSETW)+KGL-1

  ! CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G%NLOEN(IGLG),KFIELD,KFIELD)
  ! !$ACC HOST_DATA USE_DEVICE(PREEL_REAL,PREEL_COMPLEX)
  ! CALL EXECUTE_PLAN_FFTC(IPLAN_R2C,-1,PREEL_REAL(1,IOFF_REAL),PREEL_COMPLEX(1,IOFF_COMPLEX))
  ! !$ACC END HOST_DATA
END DO
IRET = CUDA_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL GSTATS(433,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(433,1)
ENDIF
CALL GSTATS(413,1)

!$ACC END DATA

!$ACC EXIT DATA DELETE(PREEL_REAL)
DEALLOCATE(PREEL_REAL)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
