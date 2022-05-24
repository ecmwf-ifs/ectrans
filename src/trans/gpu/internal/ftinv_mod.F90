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
SUBROUTINE FTINV(PREEL,KSTRIDE,KFIELDS)

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

USE TPM_DISTR       ,ONLY : D,D_NSTAGTF,D_NPTRLS, MYSETW,  MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G, G_NLOEN, G_NMEN
USE TPM_GEN         ,ONLY : NOUT, LSYNC_TRANS
USE TPM_HICFFT      ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE DEVICE_MOD
USE ISO_C_BINDING

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KSTRIDE
REAL(KIND=JPRBT), INTENT(INOUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG, IOFF, KGL, JF, JJ
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time
TYPE(C_PTR) :: IPLAN_C2R
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISIZE, IDIM2
INTEGER :: ISTAT,IDEV, IUNIT
INTEGER :: I, J

REAL(KIND=JPRBT), ALLOCATABLE  :: ZREEL2(:,:)

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
IDIM2=size(PREEL,2)
ALLOCATE(ZREEL2(ISIZE,IDIM2))

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTINV BARRIER')
ENDIF
CALL GSTATS(451,0)

#ifdef ACCGPU
!$ACC DATA CREATE(ZREEL2) PRESENT(PREEL)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(PRESENT,ALLOC:ZREEL2,PREEL)
#endif

DO KGL=IBEG,IEND,IINC
  IOFF=D_NSTAGTF(KGL)+1
  IGLG  = D_NPTRLS(MYSETW)+KGL-1
  CALL CREATE_PLAN_FFT(IPLAN_C2R,1,G%NLOEN(IGLG),2*KFIELDS,KSTRIDE)
  CALL EXECUTE_PLAN_FFT(1,G_NLOEN(IGLG),PREEL(1, IOFF),ZREEL2(1, IOFF),IPLAN_C2R)
END DO

ISTAT = DEVICE_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTINV BARRIER')
ENDIF
CALL GSTATS(451,1)

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(IGLG,JJ) DEFAULT(NONE) &
!$OMP&    SHARED(ZREEL2,PREEL,IDIM2,ISIZE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(JF,JJ) DEFAULT(NONE) &
!$ACC&    COPYIN(ISIZE,IDIM2) &
!$ACC&    PRESENT(ZREEL2,PREEL)
#endif
DO JJ=1, IDIM2
  DO JF=1,ISIZE
    PREEL(JF,JJ) = ZREEL2(JF, JJ)
  ENDDO
ENDDO
!     ------------------------------------------------------------------
!

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif
DEALLOCATE(ZREEL2)

!!$ACC UPDATE HOST(PREEL)
!!write(301,*) 'debug ftinv: ',PREEL(1,1),PREEL(1,2),PREEL(1,3)

END SUBROUTINE FTINV
END MODULE FTINV_MOD
