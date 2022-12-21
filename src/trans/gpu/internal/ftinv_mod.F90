! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
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

USE PARKIND_ECTRANS ,ONLY : JPIM, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC, D_NSTAGTF, D_NPTRLS
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN, G_NLOEN
USE TPM_GEN         ,ONLY : NOUT
USE TPM_FFT         ,ONLY : T
USE TPM_FFTH        ,ONLY : CREATE_PLAN_FFT, destroy_plan_fft, EXECUTE_PLAN_FFT
USE TPM_DIM         ,ONLY : R, R_NNOEXTZL
USE HIP_DEVICE_MOD
USE ISO_C_BINDING

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM) :: KGL
REAL(KIND=JPRBT), INTENT(INOUT), TARGET  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time
TYPE(C_PTR) :: IPLAN_C2R
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISIZE
integer :: istat,idev

REAL(KIND=JPRBT), allocatable, target  :: PREEL2(:,:)

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

#ifdef ACCGPU
!$ACC DATA &
!$ACC& PRESENT(PREEL)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:PREEL)
#endif

!stop "DEBUGGING: beginning of ftinv_mod"

#ifdef ACCGPU
!$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IOFF,IGLG,IST,ILEN,IST1) &
!$ACC&     COPYIN(IBEG,IEND,IINC,KFIELDS,MYSETW) &
!$ACC&     PRESENT(D,D%NSTAGTF,D%NPTRLS,G,G%NMEN,G%NLOEN,R,R%NNOEXTZL,PREEL)
#endif
#ifdef OMPGPU
!$OMP TARGET PARALLEL DO DEFAULT(NONE) PRIVATE(IOFF,IGLG,IST,ILEN,IST1) &
!$OMP&     SHARED(IBEG,IEND,IINC,D,MYSETW,G,R,KFIELDS,PREEL)
#endif
DO KGL=IBEG,IEND,IINC

  IOFF  = D%NSTAGTF(KGL)+1
  IGLG  = D%NPTRLS(MYSETW)+KGL-1
  IST   = 2*(G%NMEN(IGLG)+1)
  ILEN  = G%NLOEN(IGLG)+R%NNOEXTZL+2-IST
  IST1=1
  IF (G%NLOEN(IGLG)==1) IST1=0

#ifdef ACCGPU
  !$ACC LOOP COLLAPSE(2)
#endif
  !!$OMP TARGET PARALLEL DO COLLAPSE(2)
  !This OMP TARGET statement doesn't compile with nvhpc
  DO JJ=IST1,ILEN
     DO JF=1,KFIELDS
        PREEL(JF,IST+IOFF+JJ-1) = 0.0_JPRBT
     ENDDO
  ENDDO

END DO
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif

ALLOCATE(PREEL2(SIZE(PREEL,1),SIZE(PREEL,2)))
#ifdef ACCGPU
!$ACC DATA CREATE(PREEL2) PRESENT(PREEL)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:PREEL2,PREEL)
#endif

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(istat,KGL,IOFF,IGLG,IPLAN_C2R)
DO KGL=IBEG,IEND,IINC
  IOFF=D%NSTAGTF(KGL)+1
  IGLG  = D%NPTRLS(MYSETW)+KGL-1
  !IF (G%NLOEN(IGLG)>1) THEN
!call cudaProfilerStop()
     !istat=cuda_SetDevice(idev)
     CALL CREATE_PLAN_FFT(IPLAN_C2R,1,G%NLOEN(IGLG),KFIELDS)
     !!!!DEBUGGING:
     !!!!$ACC PARALLEL LOOP
     !!!DO JJ=1,10
     !!!   PREEL(JJ,ioff) = JJ
     !!!END DO
     !!!!$ACC UPDATE HOST(PREEL)
     !!!print*,"DEBUGGING: ftinv_mod: data before FFT: ",PREEL(1:10,ioff)
#ifdef ACCGPU
     ! moved this HOST_DATA USE_DEVICE into tpm_ffth.F90 to make it work on AMD GPUs (Andreas)
     !!$ACC HOST_DATA USE_DEVICE(PREEL,PREEL2)
#endif
#ifdef OMPGPU
     !!$OMP TARGET DATA USE_DEVICE_PTR(PREEL,PREEL2)
#endif
     !stop "DEBUGGING: before first execute_plan"
     CALL EXECUTE_PLAN_FFT(1,G%NLOEN(IGLG),PREEL(1, ioff),PREEL2(1, ioff),IPLAN_C2R)
     !!!CALL EXECUTE_PLAN_FFT(1,10,PREEL(1, ioff),PREEL2(1, ioff),IPLAN_C2R)
#ifdef OMPGPU
     !!$OMP END TARGET DATA
#endif
!!!!$ACC UPDATE HOST(PREEL2)
!!!print*,"DEBUGGING: ftinv_mod: data after FFT: ",PREEL2(1:10,ioff)
!!!!stop "after first execute_plan"
#ifdef ACCGPU
     !!$ACC END HOST_DATA
#endif
!call cudaProfilerStart()
  !ENDIF
END DO
!!$OMP END PARALLEL DO
ISTAT = HIP_SYNCHRONIZE()      


#ifdef OMPGPU
!$OMP TARGET
#endif
#ifdef ACCGPU
!$ACC KERNELS
#endif
PREEL(:,:) = PREEL2(:,:)
#ifdef ACCGPU
!$ACC END KERNELS
#endif
#ifdef OMPGPU
!$OMP END TARGET
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif
!     ------------------------------------------------------------------
!stop "end of ftinv_mod"
END SUBROUTINE FTINV
END MODULE FTINV_MOD
