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
SUBROUTINE FTINV(PREELIN,PREELOUT,KFIELDS)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

!        Explicit arguments :  PREELTMP  - in Fourier/grid-point array
!                              PREEL     - out Fourier/grid-point array  
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

USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_GEN         ,ONLY : NOUT
USE TPM_FFT         ,ONLY : T
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT, DESTROY_PLAN_FFT
USE TPM_DIM         ,ONLY : R
USE CUDA_DEVICE_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM) :: KGL
REAL(KIND=JPRBT), INTENT(INOUT)  :: PREELIN(:,:)
REAL(KIND=JPRBT), INTENT(OUT)  :: PREELOUT(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time
INTEGER(KIND=JPIM) :: IPLAN_C2R
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISIZE
INTEGER :: ISTAT,IDEV

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

ISIZE=size(PREELIN,1)

!$ACC DATA &
!$ACC& COPYIN(D,D%NSTAGTF,D%NPTRLS,G%NMEN,G%NLOEN,R,R%NNOEXTZL) &
!$ACC& COPYIN(D%NSTAGTF,D%NPTRLS,G%NMEN,G%NLOEN,R%NNOEXTZL) &
!$ACC& PRESENT(PREELIN)

!$ACC PARALLEL LOOP
DO KGL=IBEG,IEND,IINC

  IOFF  = D%NSTAGTF(KGL)+1
  IGLG  = D%NPTRLS(MYSETW)+KGL-1
  IST   = 2*(G%NMEN(IGLG)+1)
  ILEN  = G%NLOEN(IGLG)+R%NNOEXTZL+2-IST
  IST1=1
  IF (G%NLOEN(IGLG)==1) IST1=0

  !$ACC LOOP COLLAPSE(2)
  DO JJ=IST1,ILEN
     DO JF=1,KFIELDS
        PREELIN(JF,IST+IOFF+JJ-1) = 0.0_JPRBT
     ENDDO
  ENDDO

END DO
!$ACC END DATA

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ISTAT,KGL,IOFF,IGLG,IPLAN_C2R)
DO KGL=IBEG,IEND,IINC
  IOFF=D%NSTAGTF(KGL)+1
  IGLG  = D%NPTRLS(MYSETW)+KGL-1
  !IF (G%NLOEN(IGLG)>1) THEN
     CALL CREATE_PLAN_FFT(IPLAN_C2R,1,G%NLOEN(IGLG),KFIELDS)
     !$ACC HOST_DATA USE_DEVICE(PREELIN,PREELOUT)
     CALL EXECUTE_PLAN_FFTC(IPLAN_C2R,1,PREELIN(1, IOFF),PREELOUT(1, IOFF))
     !$ACC END HOST_DATA
  !ENDIF
END DO
!!$OMP END PARALLEL DO
ISTAT = CUDA_SYNCHRONIZE()      

!     ------------------------------------------------------------------

END SUBROUTINE FTINV
END MODULE FTINV_MOD
