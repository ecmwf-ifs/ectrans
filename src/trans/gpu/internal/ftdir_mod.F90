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

USE PARKIND_ECTRANS ,ONLY : JPIM, JPIB, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC,D_NSTAGTF,D_NPTRLS
USE TPM_TRANS       ,ONLY : ZGTF
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX
USE TPM_FFT         ,ONLY : T
USE TPM_FFTH        ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE TPM_DIM         ,ONLY : R,R_NNOEXTZL
USE HIP_DEVICE_MOD
USE ISO_C_BINDING
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

OFFSET_VAR=D_NPTRLS(MYSETW)

IMAX = G_NLOEN_MAX + 2 + R_NNOEXTZL

IDIM1=size(zgtf,1)
IDIM2=size(zgtf,2)
ALLOCATE(ZGTF2(IDIM1,IDIM2))
#ifdef ACCGPU
!$ACC ENTER DATA CREATE(ZGTF2)
!$ACC KERNELS
ZGTF2(:,:) = 0._JPRBT
!$ACC END KERNELS
#endif
!TODO: This needs to be checked by an OpenMP expert (or even an OpenMP rookie)
#ifdef OMPGPU
!$OMP TARGET DATA MAP(ALLOC:ZGTF2)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
DO I = 1, IDIM1
  DO J = 1, IDIM2
    ZGTF2(I,J) = 0._JPRBT
  END DO
END DO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif


DO KGL=IBEG,IEND,IINC
  IOFF=D_NSTAGTF(KGL)+1
  IGLG = D_NPTRLS(MYSETW)+KGL-1

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G_NLOEN(IGLG),KFIELDS)
  CALL EXECUTE_PLAN_FFT(-1,G_NLOEN(IGLG),ZGTF(1,IOFF),ZGTF2(1,IOFF),IPLAN_R2C)
END DO

ISTAT = HIP_SYNCHRONIZE()

! need a faster way for this, in place transforms ? Nils
#ifdef ACCGPU
!$ACC DATA
#endif

#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(2) PRIVATE(IGLG,JJ) DEFAULT(NONE) &
!$OMP&    SHARED(ZGTF,ZGTF2,IDIM1,IDIM2)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(JF,JJ) DEFAULT(NONE) &
!$ACC&    COPYIN(IDIM1,IDIM2) &
!$ACC&    PRESENT(ZGTF,ZGTF2)
#endif
DO JJ=1, IDIM2
  DO JF=1,IDIM1
    ZGTF(JF,JJ) = ZGTF2(JF, JJ)
  ENDDO
ENDDO

#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(3) PRIVATE(JMAX,KGL,IOFF,SCAL,IST) DEFAULT(NONE) &
!$OMP&    SHARED(IBEG,OFFSET_VAR,IEND,IINC,IMAX,KFIELDS,G_NLOEN,G_NMEN,D_NSTAGTF,ZGTF,R_NNOEXTZL)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(JMAX,KGL,IOFF,SCAL,IST) DEFAULT(NONE) &
!$ACC&    COPYIN(IBEG,IEND,IINC,OFFSET_VAR,KFIELDS,IMAX) &
!$ACC&    PRESENT(ZGTF,G_NLOEN,G_NMEN,D_NSTAGTF,R_NNOEXTZL)
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
           ZGTF(JF,IOFF+JJ-1)= SCAL * ZGTF(JF, IOFF+JJ-1)
         END IF

         ! case JJ>0
         IF( JJ .LE. (JMAX+R_NNOEXTZL+2-IST)) ZGTF(JF,IST+IOFF+JJ-1) = 0.0_JPRBT
         ! case JJ=0
         IF (G_NLOEN(IGLG)==1) ZGTF(JF,IST+IOFF-1) = 0.0_JPRBT
      ENDDO
   ENDDO
ENDDO

#ifdef ACCGPU
!$ACC END DATA
#endif

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif

#ifdef ACCGPU
!$ACC EXIT DATA DELETE(ZGTF2)
#endif
DEALLOCATE(ZGTF2)
!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
