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

MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(ZGTF,KSTRIDE,KFIELD)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  ZGTF    - Fourier/grid-point array
!        --------------------  KSTRIDE - stride of PREEL
!                              KFIELD   - number of fields

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

USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC,D_NSTAGTF,D_NPTRLS, D_NSTAGT0B, D_NPNTGTB0, D_NPROCM
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NLOEN
USE TPM_FFT         ,ONLY : T
USE TPM_HICFFT      ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE DEVICE_MOD
USE ISO_C_BINDING
use openacc
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KSTRIDE, KFIELD
REAL(KIND=JPRBT), INTENT(INOUT) :: ZGTF(:,:)
INTEGER(KIND=JPIM)  :: KGL

INTEGER(KIND=JPIM) :: IGLG,JM,JF,IST1,IPROC,ISTA
INTEGER(KIND=JPIM) :: IOFF,D_NDGL_FS
TYPE(C_PTR)        :: IPLAN_R2C

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,IRET,D_NDGL_FS
INTEGER(KIND=JPIM) :: ISIZE, II, IMAX, IDIM1, IDIM2, I, J
REAL(KIND=JPRBT), ALLOCATABLE :: ZGTF2(:,:)


!     ------------------------------------------------------------------

D_NDGL_FS=D%NDGL_FS

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D_NDGL_FS
  IINC=1
ELSE
  IBEG=D_NDGL_FS
  IEND=1
  IINC=-1
ENDIF

D_NDGL_FS = D%NDGL_FS

IDIM1=SIZE(ZGTF,1)
IDIM2=SIZE(ZGTF,2)
ALLOCATE(ZGTF2(IDIM1,IDIM2))

#ifdef ACCGPU
!$ACC DATA PRESENT(ZGTF, &
!$ACC&             D_NSTAGTF,D_NSTAGT0B,D_NPTRLS,D_NPROCM,D_NPNTGTB0,G_NMEN, &
!$ACC&             G_NLOEN) &
!$ACC&     CREATE(ZGTF2)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:ZGTF2)
#endif

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTDIR BARRIER')
ENDIF
CALL GSTATS(450,0)

DO KGL=IBEG,IEND,IINC

  IOFF=D_NSTAGTF(KGL)+1
  IGLG = D_NPTRLS(MYSETW)+KGL-1

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G_NLOEN(IGLG),KFIELD,KFIELD)
  CALL EXECUTE_PLAN_FFT(-1,G_NLOEN(IGLG),ZGTF(1,IOFF),ZGTF2(1,IOFF),IPLAN_R2C)
ENDDO

IRET = DEVICE_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTDIR BARRIER')
ENDIF
CALL GSTATS(450,1)

#ifdef ACCGPU
!$ACC KERNELS DEFAULT(NONE)
#endif
ZGTF(:,:) = ZGTF2(:,:)
#ifdef ACCGPU
!$ACC END KERNELS
#endif


#ifdef ACCGPU
!$ACC END DATA
#endif

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif


DEALLOCATE(ZGTF2)
!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
