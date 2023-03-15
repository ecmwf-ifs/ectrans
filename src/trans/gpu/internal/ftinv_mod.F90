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

USE TPM_DISTR       ,ONLY : D,D_NSTAGTF,D_NPTRLS, MYSETW,  MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G, G_NLOEN, G_NMEN
USE TPM_GEN         ,ONLY : NOUT
USE TPM_FFT         ,ONLY : T
USE TPM_FFTH        ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE TPM_DIM         ,ONLY : R, R_NNOEXTZL
USE HIP_DEVICE_MOD
USE ISO_C_BINDING

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM) :: KGL
REAL(KIND=JPRBT), INTENT(INOUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time
TYPE(C_PTR) :: IPLAN_C2R
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISIZE, IDIM2
integer :: istat,idev, iunit

REAL(KIND=JPRBT), allocatable  :: ZREEL2(:,:)

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
!WARNING: following ALLOC statement should be PRESENT,ALLOC but cause issues with AMD compiler!
!$OMP TARGET DATA MAP(ALLOC:PREEL)
#endif

#ifdef ACCGPU
!$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IOFF,IGLG,IST,ILEN,IST1) &
!$ACC&     COPYIN(IBEG,IEND,IINC,KFIELDS,MYSETW) &
!$ACC&     PRESENT(D_NSTAGTF,D_NPTRLS,G_NMEN,G_NLOEN,R,R_NNOEXTZL,PREEL)
#endif
#ifdef OMPGPU
!WARNING: in OpenACC we have nested ACC loops. Can't get this to work in OpenMP with AMD compiler
!!$OMP TARGET PARALLEL DO DEFAULT(NONE) PRIVATE(IOFF,IGLG,IST,ILEN,IST1) &
!!$OMP&     SHARED(IBEG,IEND,IINC,MYSETW,KFIELDS,PREEL,D_NSTAGTF,D_NPTRLS,G_NMEN,G_NLOEN,R,R_NNOEXTZL)
#endif
DO KGL=IBEG,IEND,IINC

  IOFF  = D_NSTAGTF(KGL)+1
  IGLG  = D_NPTRLS(MYSETW)+KGL-1
  IST   = 2*(G_NMEN(IGLG)+1)
  ILEN  = G_NLOEN(IGLG)+R_NNOEXTZL+2-IST
  IST1=1
  IF (G_NLOEN(IGLG)==1) IST1=0

#ifdef ACCGPU
  !$ACC LOOP COLLAPSE(2)
#endif
#ifdef OMPGPU
  !$OMP TARGET PARALLEL DO DEFAULT(NONE) PRIVATE(JJ,JF) SHARED(ILEN,KFIELDS,PREEL,IST1,IST,IOFF)
#endif
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

IDIM2=size(PREEL,2)
ALLOCATE(ZREEL2(ISIZE,IDIM2))
#ifdef ACCGPU
!$ACC ENTER DATA CREATE(ZREEL2)
!$ACC KERNELS
ZREEL2(:,:) = 0._JPRBT
!$ACC END KERNELS
#endif
!TODO: This needs to be checked by an OpenMP expert (or even an OpenMP rookie)
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:ZREEL2,PREEL)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
ZREEL2(:,:) = 0._JPRBT
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif

DO KGL=IBEG,IEND,IINC
  IOFF=D_NSTAGTF(KGL)+1
  IGLG  = D_NPTRLS(MYSETW)+KGL-1
  CALL CREATE_PLAN_FFT(IPLAN_C2R,1,G_NLOEN(IGLG),KFIELDS)
  CALL EXECUTE_PLAN_FFT(1,G_NLOEN(IGLG),PREEL(1, ioff),ZREEL2(1, ioff),IPLAN_C2R)
END DO
ISTAT = HIP_SYNCHRONIZE()

#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(2) PRIVATE(IGLG,JJ) DEFAULT(NONE) &
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

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC EXIT DATA DELETE(ZREEL2)
#endif
DEALLOCATE(ZREEL2)

END SUBROUTINE FTINV
END MODULE FTINV_MOD
