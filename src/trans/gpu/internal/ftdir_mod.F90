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
SUBROUTINE FTDIR(KFIELDS)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

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
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS ,ONLY : JPIM, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC,D_NSTAGTF,D_NPTRLS
USE TPM_TRANS       ,ONLY : ZGTF
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX
USE TPM_FFT         ,ONLY : T
USE TPM_HICFFT      ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE TPM_DIM         ,ONLY : R,R_NNOEXTZL
USE DEVICE_MOD
USE ISO_C_BINDING
use openacc
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS
INTEGER(KIND=JPIM)  :: KGL
!!!REAL(KIND=JPRBT), INTENT(INOUT) :: PREEL(KFIELDS,D%NLENGTF)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
TYPE(C_PTR) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: JMAX
REAL(KIND=JPRBT)    :: SCAL
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISCAL
INTEGER(KIND=JPIM) :: OFFSET_VAR, IUNIT, ISIZE, II, IMAX, IDIM1, IDIM2, I, J
integer :: istat, idev
REAL(KIND=JPRBT), ALLOCATABLE :: ZGTF2(:,:)

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

IDIM1=SIZE(ZGTF,1)
IDIM2=SIZE(ZGTF,2)
ALLOCATE(ZGTF2(IDIM1,IDIM2))

#ifdef ACCGPU
!$ACC DATA PRESENT(ZGTF,D_NSTAGTF,D_NPTRLS,G_NMEN, &
!$ACC&             G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX,R_NNOEXTZL)
!$ACC DATA CREATE(ZGTF2)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:ZGTF2)
#endif


!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(KGL,IOFF,IGLG,IPLAN_R2C,istat)
DO KGL=IBEG,IEND,IINC

  IOFF=D_NSTAGTF(KGL)+1
  IGLG = D_NPTRLS(MYSETW)+KGL-1

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G_NLOEN(IGLG),KFIELDS)
  CALL EXECUTE_PLAN_FFT(-1,G_NLOEN(IGLG),ZGTF(1,IOFF),ZGTF2(1,IOFF),IPLAN_R2C)
ENDDO
!!$OMP END PARALLEL DO

ISTAT = DEVICE_SYNCHRONIZE()

OFFSET_VAR=D_NPTRLS(MYSETW)
#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(JMAX,KGL,IOFF,SCAL,IST) DEFAULT(NONE) &
!$OMP&    SHARED(IBEG,OFFSET_VAR,IEND,IINC,IMAX,KFIELDS,G_NLOEN,G_NMEN,D_NSTAGTF,ZGTF,R_NNOEXTZL)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(JMAX,JJ,KGL,IOFF,SCAL,IST) FIRSTPRIVATE(IBEG,IEND,IINC,OFFSET_VAR,KFIELDS)
#endif
DO IGLG=IBEG+OFFSET_VAR-1,IEND+OFFSET_VAR-1,IINC
  DO JF=1,KFIELDS
    JMAX = G_NLOEN(IGLG)
    SCAL = 1._JPRBT/REAL(JMAX,JPRBT)
    IST  = 2*(G_NMEN(IGLG)+1)
    KGL=IGLG-OFFSET_VAR+1
    IOFF=D_NSTAGTF(KGL)+1

    !$ACC LOOP SEQ
    DO JJ=1, JMAX
      ZGTF(JF,IOFF+JJ-1)= SCAL * ZGTF2(JF, IOFF+JJ-1)
    ENDDO

    !! WHAT'S GOING ON HERE? TRUNCATING?
    IF (JMAX== 1) ZGTF(JF,IST+IOFF-1) = 0.0_JPRBT
    !$ACC LOOP SEQ
    DO JJ=1,JMAX+R_NNOEXTZL+3-IST
      ZGTF(JF,IST+IOFF+JJ-1) = 0.0_JPRBT
    ENDDO
  ENDDO
ENDDO

#ifdef ACCGPU
!$ACC END DATA
!$ACC END DATA
#endif

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif


DEALLOCATE(ZGTF2)
!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
