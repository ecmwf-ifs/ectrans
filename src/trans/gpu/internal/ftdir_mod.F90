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
SUBROUTINE FTDIR(STRIDE,KF_FS)


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

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, D_NSTAGTF,D_NPTRLS, D_NSTAGT0B, D_NPNTGTB0, D_NPROCM
USE TPM_TRANS       ,ONLY : ZGTF, FOUBUF_IN
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NLOEN
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE CUDA_DEVICE_MOD
USE MPL_MODULE      ,ONLY : MPL_BARRIER
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: STRIDE,KF_FS
INTEGER(KIND=JPIM)  :: KGL

INTEGER(KIND=JPIM) :: IGLG,JM,JF,IST1, IPROC, ISTA
INTEGER(KIND=JPIM) :: IOFF
INTEGER(KIND=JPIM) :: IPLAN_R2C
REAL(KIND=JPRBT)    :: SCAL

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
INTEGER(KIND=JPIM) :: OFFSET_VAR
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
!$ACC& PRESENT(ZGTF,FOUBUF_IN, &
!$ACC&         D,D_NSTAGTF,G_NMEN,G_NLOEN,D_NSTAGT0B,D_NPNTGTB0,D_NPROCM) &
!$ACC& CREATE(ZGTF2)

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

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='FTDIR BARRIER')
ENDIF
CALL GSTATS(450,1)

! scale results and move into next transformation buffer

OFFSET_VAR=D_NPTRLS(MYSETW)
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IGLG,IOFF,SCAL,JM) DEFAULT(NONE)
DO KGL=1,D%NDGL_FS
  DO JF=1,KF_FS
    IGLG = OFFSET_VAR+KGL-1
    IOFF = D_NSTAGTF(KGL)+1

    SCAL = 1._JPRBT/REAL(G_NLOEN(IGLG),JPRBT)
    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      IPROC = D_NPROCM(JM)
      ISTA  = (D_NSTAGT0B(IPROC)+D_NPNTGTB0(JM,KGL))*2*KF_FS

      FOUBUF_IN(ISTA+2*JF-1) = SCAL * ZGTF2(2*JF-1, 2*JM+IOFF)
      FOUBUF_IN(ISTA+2*JF  ) = SCAL * ZGTF2(2*JF  , 2*JM+IOFF)
     ENDDO
   ENDDO
ENDDO

!$ACC END DATA

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
