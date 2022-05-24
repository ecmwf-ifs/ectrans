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
SUBROUTINE FTDIR(ZGTF,STRIDE,KF_FS)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KSTIRDE - stride of PREEL
!                              KF_FS   - number of fields

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
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: STRIDE,KF_FS
REAL(KIND=JPRBT), INTENT(INOUT) :: ZGTF(:,:)
INTEGER(KIND=JPIM)  :: KGL

INTEGER(KIND=JPIM) :: IGLG,JM,JF,IST1, IPROC, ISTA
INTEGER(KIND=JPIM) :: IOFF
INTEGER(KIND=JPIM) :: IPLAN_R2C

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC, IRET
real(kind=jprbt), allocatable :: zgtf2(:,:)

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


allocate(zgtf2(size(zgtf,1),size(zgtf,2)))
!$ACC DATA CREATE(ZGTF2) PRESENT(ZGTF)

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTDIR BARRIER')
ENDIF
CALL GSTATS(450,0)

DO KGL=IBEG,IEND,IINC

  IOFF=D%NSTAGTF(KGL)+1
  IGLG = D%NPTRLS(MYSETW)+KGL-1

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G%NLOEN(IGLG),2*KF_FS,STRIDE)
  !$ACC host_data use_device(ZGTF,ZGTF2)
  CALL EXECUTE_PLAN_FFTC(IPLAN_R2C,-1,ZGTF(1,IOFF),ZGTF2(1,IOFF))
  !$ACC end host_data
END DO

IRET = CUDA_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTDIR BARRIER')
ENDIF
CALL GSTATS(450,1)

!$ACC KERNELS DEFAULT(NONE)
ZGTF(:,:) = ZGTF2(:,:)
!$ACC END KERNELS

!$ACC END DATA

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
