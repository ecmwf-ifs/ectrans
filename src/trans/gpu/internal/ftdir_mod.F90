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
SUBROUTINE FTDIR(PREEL_REAL,PREEL_COMPLEX,KFIELD)

!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELD  - number of fields

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
USE TPM_HICFFT      ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT, EXECUTE_DIR_FFT
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE TPM_STATS       ,ONLY : GSTATS => GSTATS_NVTX
USE DEVICE_MOD
USE ISO_C_BINDING
use openacc
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_REAL(:)
REAL(KIND=JPRBT), INTENT(OUT), POINTER :: PREEL_COMPLEX(:)

INTEGER(KIND=JPIM) :: IGLG,IOFF_REAL,IOFF_COMPLEX,KGL,IRET,D_NDGL_FS
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
TYPE(C_PTR)        :: IPLAN_R2C
REAL(KIND=JPRBT), POINTER  :: PREEL2(:,:), TMP(:,:)

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

ALLOCATE(PREEL_COMPLEX(KFIELD*D%NLENGTF))
 
#ifdef ACCGPU
!$ACC ENTER DATA CREATE(PREEL_COMPLEX)
!$ACC DATA PRESENT(PREEL_REAL, PREEL_COMPLEX, &
!$ACC&             D_NSTAGTF,D_NSTAGT0B,D_NPTRLS,D_NPROCM,D_NPNTGTB0,G_NMEN,G_NLOEN)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:PREEL_COMPLEX)
#endif

IF (LSYNC_TRANS) THEN
  CALL GSTATS(430,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(430,1)
ENDIF
CALL GSTATS(413,0)
CALL EXECUTE_DIR_FFT(PREEL_REAL(:),PREEL_COMPLEX(:),KFIELD, &
    & LOENS=G%NLOEN(D%NPTRLS(MYSETW):D%NPTRLS(MYSETW)+D%NDGL_FS-1), &
    & OFFSETS=D%NSTAGTF(1:D%NDGL_FS))

!!DO KGL=IBEG,IEND,IINC
!!
!!  ! NSTAGTF gives us space for NLOEN+3 elements
!!  ! In reality, at this point we need space for at most NLOEN+2 elements
!!  ! (in case NLOEN is even, otherwise NLOEN+1, due to the R2C definition)
!!  IOFF_REAL=D%NSTAGTF(KGL)+1
!!  IOFF_COMPLEX=D%NSTAGTF(KGL)/2+1
!!  IGLG = D_NPTRLS(MYSETW)+KGL-1
!!
!!  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G_NLOEN(IGLG),KFIELD,KFIELD)
!!  CALL EXECUTE_PLAN_FFT(-1,G_NLOEN(IGLG),PREEL_REAL(1,IOFF_REAL),PREEL_COMPLEX(1,IOFF_COMPLEX),IPLAN_R2C)
!!ENDDO

IRET = DEVICE_SYNCHRONIZE()

IF (LSYNC_TRANS) THEN
  CALL GSTATS(433,0)
  CALL MPL_BARRIER(CDSTRING='')
  CALL GSTATS(433,1)
ENDIF
CALL GSTATS(413,1)


#ifdef ACCGPU
!$ACC END DATA
#endif

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif


! and deallocate the local pointer
!$ACC EXIT DATA DELETE(PREEL_REAL)
DEALLOCATE(PREEL_REAL)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
