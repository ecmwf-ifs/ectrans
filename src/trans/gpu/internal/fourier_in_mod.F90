! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_IN_MOD
CONTAINS
SUBROUTINE FOURIER_IN(FOUBUF,PREEL,KFIELDS)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM, JPRB, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, D_NSTAGTF,D_NSTAGT0B,D_NPNTGTB0,D_NPROCM,D_NPTRLS
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NMEN_MAX
use tpm_gen, only: nout
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS

INTEGER(KIND=JPIM) :: KGL

REAL(KIND=JPRBT), INTENT(OUT) :: PREEL(:,:)
REAL(KIND=JPRB), ALLOCATABLE, INTENT(IN) :: FOUBUF(:)

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,IR,II,ISTA,OFFSET_VAR,IOFF
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,iimax1,iimax2,iimax3,iunit

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

!$ACC DATA PRESENT(G_NMEN,D_NPROCM,D_NSTAGT0B,D_NPNTGTB0,FOUBUF,PREEL,D_NSTAGTF) ASYNC(1)

! TODO: We don't need to zero out the full array here but we need to zero out because implicit
! truncation happens. We cannot rely on previous iterations that they had the same configuration.

!$ACC KERNELS ASYNC(1)
PREEL(:,:) = 0
!$ACC END KERNELS

OFFSET_VAR=D_NPTRLS(MYSETW)
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IGLG,IOFF,JM,IPROC,ISTA) DEFAULT(NONE) ASYNC(1)
DO KGL=IBEG,IEND,IINC
  DO JF=1,KFIELDS
    IGLG = OFFSET_VAR+KGL-1
    IOFF = D_NSTAGTF(KGL)+1

    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      IPROC = D_NPROCM(JM)
      ISTA  = (D_NSTAGT0B(IPROC)+D_NPNTGTB0(JM,KGL))*2*KFIELDS

      PREEL(2*JF-1,2*JM+IOFF) = FOUBUF(ISTA+2*JF-1)
      PREEL(2*JF,  2*JM+IOFF) = FOUBUF(ISTA+2*JF  )
    ENDDO
  ENDDO
ENDDO
!$ACC END DATA

!$ACC WAIT(1)

!$ACC EXIT DATA DELETE(FOUBUF)
DEALLOCATE(FOUBUF)

END SUBROUTINE FOURIER_IN
END MODULE FOURIER_IN_MOD

