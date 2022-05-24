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
SUBROUTINE FTINV(PREEL,KFIELDS)

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
use tpm_gen, only: nout
USE TPM_FFT         ,ONLY : T
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT, destroy_plan_fft
USE TPM_DIM         ,ONLY : R
USE CUDA_DEVICE_MOD
USE MPL_MODULE      ,ONLY : MPL_BARRIER

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM) :: KGL
REAL(KIND=JPRBT), INTENT(INOUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time
INTEGER(KIND=JPIM) :: IPLAN_C2R
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISIZE
integer :: istat,idev

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

ISIZE=size(PREEL,1)

!$ACC DATA &
!$ACC& PRESENT(PREEL)

!$ACC PARALLEL LOOP DEFAULT(NONE)
DO KGL=IBEG,IEND,IINC

  IOFF  = D%NSTAGTF(KGL)+1
  IGLG  = D%NPTRLS(MYSETW)+KGL-1
  IST   = 2*(G%NMEN(IGLG)+1)
  ILEN  = G%NLOEN(IGLG)+R%NNOEXTZL+2-IST
  IST1=1
  IF (G%NLOEN(IGLG)==1) IST1=0

  !$ACC loop collapse(2)
  DO JJ=IST1,ILEN
     DO JF=1,KFIELDS
        PREEL(JF,IST+IOFF+JJ-1) = 0.0_JPRBT
     ENDDO
  ENDDO

END DO
!$ACC end data

allocate(preel2(size(preel,1),size(preel,2)))
!$acc data create(preel2) present(preel)

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTINV BARRIER')
ENDIF
CALL GSTATS(451,0)

!istat = cuda_GetDevice(idev)
!istat = cuda_Synchronize()      
!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(istat,KGL,IOFF,IGLG,IPLAN_C2R)
DO KGL=IBEG,IEND,IINC
  IOFF=D%NSTAGTF(KGL)+1
  IGLG  = D%NPTRLS(MYSETW)+KGL-1
  !IF (G%NLOEN(IGLG)>1) THEN
!call cudaProfilerStop()
     !istat=cuda_SetDevice(idev)
     CALL CREATE_PLAN_FFT(IPLAN_C2R,1,G%NLOEN(IGLG),KFIELDS,KFIELDS)
     !$ACC host_data use_device(PREEL,PREEL2)
     CALL EXECUTE_PLAN_FFTC(IPLAN_C2R,1,PREEL(1, ioff),PREEL2(1, ioff))
     !$ACC end host_data
!call cudaProfilerStart()
  !ENDIF
END DO
!!$OMP END PARALLEL DO
istat = cuda_Synchronize()      

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
