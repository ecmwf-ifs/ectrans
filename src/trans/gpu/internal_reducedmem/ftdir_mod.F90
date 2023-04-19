! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
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

USE PARKIND_ECTRANS ,ONLY : JPIM, JPIB, JPRB, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC,D_NSTAGTF,D_NPTRLS
USE TPM_TRANS       ,ONLY : ZGTF, ZGTFTMP
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX
USE TPM_FFT         ,ONLY : T
USE TPM_FFTH        ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
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
REAL(KIND=JPRBT)   :: SCAL
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISCAL
INTEGER(KIND=JPIM) :: OFFSET_VAR, IUNIT, ISIZE, II, IMAX
integer :: istat, idev

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

OFFSET_VAR=D_NPTRLS(MYSETW)

IMAX = G_NLOEN_MAX + 2 + R_NNOEXTZL

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(KGL,IOFF,IGLG,IPLAN_R2C,istat)
DO KGL=IBEG,IEND,IINC

  IOFF=D%NSTAGTF(KGL)+1
  IGLG = D%NPTRLS(MYSETW)+KGL-1
  !ILEN = G_NLOEN(IGLG)+R_NNOEXTZL+3-IST
  !IRLEN=G_NLOEN(IGLG)+R_NNOEXTZL
  !ICLEN=(IRLEN/2+1)*2

  !istat = cuda_SetDevice(idev)
  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G%NLOEN(IGLG),KFIELDS)
#ifdef ACCGPU
  !$ACC HOST_DATA USE_DEVICE(ZGTF, ZGTFTMP)
#endif
#ifdef OMPGPU
  !$OMP TARGET DATA USE_DEVICE_PTR(ZGTF, ZGTFTMP)
#endif
  CALL EXECUTE_PLAN_FFT(-1,G%NLOEN(IGLG),C_LOC(ZGTF(1,IOFF)),C_LOC(ZGTFTMP(1,IOFF)),IPLAN_R2C)
#ifdef OMPGPU
  !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
  !$ACC END HOST_DATA
#endif
END DO
!!$OMP END PARALLEL DO

ISTAT = CUDA_SYNCHRONIZE()

#ifdef ACCGPU
!$ACC DATA &
!$ACC& COPY(D,D_NSTAGTF,D_NPTRLS,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX,R_NNOEXTZL) &
!$ACC& PRESENT(ZGTFTMP)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA &
!$OMP& MAP(TO:D,D_NSTAGTF,D_NPTRLS,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX,R_NNOEXTZL) &
!$OMP& MAP(ALLOC:ZGTFTMP)
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(JMAX,KGL,IOFF,SCAL,IST)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(JMAX,KGL,IOFF,SCAL,IST)
#endif
DO IGLG=IBEG+OFFSET_VAR-1,IEND+OFFSET_VAR-1,IINC
   DO JJ=1, IMAX
      DO JF=1,KFIELDS
         JMAX = G_NLOEN(IGLG)
         IST  = 2*(G_NMEN(IGLG)+1)
         IF (JJ .LE. JMAX) THEN
           KGL=IGLG-OFFSET_VAR+1
           IOFF=D_NSTAGTF(KGL)+1
           SCAL = 1._JPRBT/REAL(G_NLOEN(IGLG),JPRBT)
           ZGTFTMP(JF,IOFF+JJ-1)= SCAL * ZGTFTMP(JF, IOFF+JJ-1)
         END IF

         ! case JJ>0
         IF( JJ .le. (JMAX+R_NNOEXTZL+2-IST)) ZGTFTMP(JF,IST+IOFF+JJ-1) = 0.0_JPRBT
         ! case JJ=0
         IF (G_NLOEN(IGLG)==1) ZGTFTMP(JF,IST+IOFF-1) = 0.0_JPRBT
      ENDDO
   ENDDO
ENDDO

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
