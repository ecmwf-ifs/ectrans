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
SUBROUTINE FTDIR(PREEL,KFIELD)

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
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE CUDA_DEVICE_MOD
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IOFF,KGL,IRET
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
REAL(KIND=JPRBT), POINTER  :: PREEL2(:,:), TMP(:,:)

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


ALLOCATE(PREEL2(SIZE(PREEL,1),SIZE(PREEL,2)))
!$ACC ENTER DATA CREATE(PREEL2)

!$ACC DATA PRESENT(PREEL,PREEL2)

IF (LSYNC_TRANS) THEN
  CALL GSTATS(430,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(430,1)
ENDIF
CALL GSTATS(413,0)

DO KGL=IBEG,IEND,IINC

  ! NSTAGTF gives us space for NLOEN+3 elements
  ! In reality, at this point we need space for at most NLOEN+2 elements
  ! (in case NLOEN is even, otherwise NLOEN+1, due to the R2C definition)
  IOFF=D%NSTAGTF(KGL)+1
  IGLG = D%NPTRLS(MYSETW)+KGL-1

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G%NLOEN(IGLG),KFIELD,KFIELD*2)
  !$ACC HOST_DATA USE_DEVICE(PREEL,PREEL2)
  CALL EXECUTE_PLAN_FFTC(IPLAN_R2C,-1,PREEL(1,IOFF),PREEL2(1,IOFF))
  !$ACC END HOST_DATA
END DO
IRET = CUDA_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL GSTATS(433,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(433,1)
ENDIF
CALL GSTATS(413,1)

!$ACC END DATA

! Swap pointers
TMP => PREEL
PREEL => PREEL2
PREEL2 => TMP

! and deallocate the local pointer
!$ACC EXIT DATA DELETE(PREEL2)
DEALLOCATE(PREEL2)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
