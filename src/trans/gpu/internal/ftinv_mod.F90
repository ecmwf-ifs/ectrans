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
SUBROUTINE FTINV(PREEL,STRIDE,KFIELDS)

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
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
!     ------------------------------------------------------------------

USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE PARKIND_ECTRANS  ,ONLY : JPIM, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE CUDA_DEVICE_MOD
USE MPL_MODULE      ,ONLY : MPL_BARRIER

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS, STRIDE
REAL(KIND=JPRBT), INTENT(INOUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IOFF,KGL
INTEGER(KIND=JPIM) :: IPLAN_C2R
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

REAL(KIND=JPRBT), allocatable  :: PREEL2(:,:)

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

allocate(preel2(size(preel,1),size(preel,2)))
!$acc data create(preel2) present(preel)

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTINV BARRIER')
ENDIF
CALL GSTATS(451,0)

DO KGL=IBEG,IEND,IINC

  IOFF=D%NSTAGTF(KGL)+1
  IGLG  = D%NPTRLS(MYSETW)+KGL-1

  CALL CREATE_PLAN_FFT(IPLAN_C2R,1,G%NLOEN(IGLG),2*KFIELDS,STRIDE)
  !$ACC HOST_DATA USE_DEVICE(PREEL,PREEL2)
  CALL EXECUTE_PLAN_FFTC(IPLAN_C2R,1,PREEL(1,IOFF),PREEL2(1,IOFF))
  !$ACC END HOST_DATA
END DO

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTINV BARRIER')
ENDIF
CALL GSTATS(451,1)

!$acc kernels
preel(:,:) = preel2(:,:)
!$acc end kernels
!$acc end data
!     ------------------------------------------------------------------

END SUBROUTINE FTINV
END MODULE FTINV_MOD
