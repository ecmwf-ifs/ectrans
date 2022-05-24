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

REAL(KIND=JPRBT), ALLOCATABLE :: PREEL_REAL_TMP(:)
INTEGER(KIND=JPIM) :: IGLG,IOFF_REAL,IOFF_COMPLEX,KGL,IRET
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
INTEGER(KIND=JPIM) :: KF, KI

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


ALLOCATE(PREEL_REAL_TMP(KFIELD*D%NLENGTF))
ALLOCATE(PREEL_COMPLEX(KFIELD*D%NLENGTF))
!$ACC ENTER DATA CREATE(PREEL_COMPLEX,PREEL_REAL_TMP)

!$ACC DATA PRESENT(PREEL_REAL,PREEL_COMPLEX,PREEL_REAL_TMP)

IF (LSYNC_TRANS) THEN
  CALL GSTATS(430,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(430,1)
ENDIF
CALL GSTATS(413,0)
!$ACC KERNELS
PREEL_REAL_TMP(:) = -1
!$ACC END KERNELS
!$ACC PARALLEL LOOP COLLAPSE(2)
DO KGL=1,D%NDGL_FS
  DO KF=1,KFIELD
    !$ACC LOOP SEQ
    DO KI=1,D%NSTAGTF(KGL+1)-D%NSTAGTF(KGL)
      PREEL_REAL_TMP(KFIELD*D%NSTAGTF(KGL)+(KF-1)*(D%NSTAGTF(KGL+1)-D%NSTAGTF(KGL))+KI) = &
          & PREEL_REAL(KFIELD*D%NSTAGTF(KGL)+(KI-1)*KFIELD+KF)
    ENDDO
  ENDDO
ENDDO
CALL EXECUTE_DIR_FFT(PREEL_REAL_TMP(:),PREEL_COMPLEX(:),KFIELD, &
    & LOENS=G%NLOEN(D%NPTRLS(MYSETW):D%NPTRLS(MYSETW)+D%NDGL_FS-1), &
    & OFFSETS=D%NSTAGTF(1:D%NDGL_FS+1))

IF (LSYNC_TRANS) THEN
  CALL GSTATS(433,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(433,1)
ENDIF
CALL GSTATS(413,1)

!$ACC END DATA

!$ACC EXIT DATA DELETE(PREEL_REAL,PREEL_REAL_TMP)
DEALLOCATE(PREEL_REAL)
DEALLOCATE(PREEL_REAL_TMP)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
