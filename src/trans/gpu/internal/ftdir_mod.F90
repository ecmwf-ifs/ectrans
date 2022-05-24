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

USE PARKIND_ECTRANS  ,ONLY : JPIM, JPIB, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC,D_NSTAGTF,D_NPTRLS
USE TPM_TRANS       ,ONLY : ZGTF
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX
USE TPM_FFT         ,ONLY : T
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE TPM_DIM         ,ONLY : R,R_NNOEXTZL
USE CUDA_DEVICE_MOD
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS
INTEGER(KIND=JPIM)  :: KGL
!!!REAL(KIND=JPRBT), INTENT(INOUT) :: PREEL(KFIELDS,D%NLENGTF)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: JMAX
REAL(KIND=JPRBT)    :: SCAL
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISCAL
INTEGER(KIND=JPIM) :: OFFSET_VAR, IUNIT, ISIZE, II, IMAX
integer :: istat, idev
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
!$ACC DATA &
!$ACC& PRESENT(ZGTF,D,D_NSTAGTF,D_NPTRLS,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX,R_NNOEXTZL)

!$ACC DATA CREATE(ZGTF2)

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(KGL,IOFF,IGLG,IPLAN_R2C,istat)
DO KGL=IBEG,IEND,IINC

  IOFF=D%NSTAGTF(KGL)+1
  IGLG = D%NPTRLS(MYSETW)+KGL-1
  !ILEN = G_NLOEN(IGLG)+R_NNOEXTZL+3-IST
  !IRLEN=G_NLOEN(IGLG)+R_NNOEXTZL
  !ICLEN=(IRLEN/2+1)*2

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G%NLOEN(IGLG),KFIELDS)
  !$ACC host_data use_device(ZGTF,ZGTF2)
  CALL EXECUTE_PLAN_FFTC(IPLAN_R2C,-1,ZGTF(1,IOFF),ZGTF2(1,IOFF))
  !$ACC end host_data
END DO
!!$OMP END PARALLEL DO

istat = cuda_Synchronize()

OFFSET_VAR=D_NPTRLS(MYSETW)
!$ACC parallel loop collapse(2) private(JMAX,JJ,KGL,IOFF,SCAL,IST) DEFAULT(NONE)
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
      DO JJ=1,JMAX+R%NNOEXTZL+3-IST
        ZGTF(JF,IST+IOFF+JJ-1) = 0.0_JPRBT
      ENDDO
   ENDDO
ENDDO
!$acc end data

!$ACC end data

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
