! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_OUT_MOD
CONTAINS
SUBROUTINE FOURIER_OUT(PREEL,FOUBUF_IN,KFIELDS)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

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

USE PARKIND_ECTRANS  ,ONLY : JPIM,JPRBT
USE TPM_DISTR       ,ONLY : D,MYSETW,MYPROC, NPROC, D_NSTAGTF, D_NSTAGT0B,D_NPNTGTB0,D_NPROCM,D_NPTRLS
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NLOEN
!

IMPLICIT NONE

REAL(KIND=JPRBT), INTENT(IN) :: PREEL(:,:)
REAL(KIND=JPRBT), ALLOCATABLE, INTENT(OUT) :: FOUBUF_IN(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,ISTA,OFFSET_VAR,IOFF,KGL
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

REAL(KIND=JPRBT)    :: SCAL

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

ALLOCATE(FOUBUF_IN(D%NLENGT0B*KFIELDS*2))
!$ACC ENTER DATA CREATE(FOUBUF_IN)

!$ACC DATA PRESENT(G_NMEN,D_NPROCM,D_NSTAGT0B,D_NPNTGTB0,FOUBUF_IN,PREEL,D_NSTAGTF,G_NLOEN) ASYNC(1)

! scale results and move into next transformation buffer

OFFSET_VAR=D_NPTRLS(MYSETW)
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IGLG,IOFF,JM,IPROC,ISTA,SCAL) DEFAULT(NONE) ASYNC(1)
DO KGL=IBEG,IEND,IINC
  DO JF=1,KFIELDS
    IGLG = OFFSET_VAR+KGL-1
    IOFF = D_NSTAGTF(KGL)+1

    SCAL = 1._JPRBT/REAL(G_NLOEN(IGLG),JPRBT)
    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      IPROC = D_NPROCM(JM)
      ISTA  = (D_NSTAGT0B(IPROC)+D_NPNTGTB0(JM,KGL))*KFIELDS*2

      ! This is not contiguous in PREEL due to the memory layout.
      FOUBUF_IN(ISTA+2*JF-1) = SCAL * PREEL(2*JF-1, 2*JM+IOFF)
      FOUBUF_IN(ISTA+2*JF  ) = SCAL * PREEL(2*JF  , 2*JM+IOFF)
    ENDDO
  ENDDO
ENDDO
!$ACC END DATA

!$ACC WAIT(1)

END SUBROUTINE FOURIER_OUT
END MODULE FOURIER_OUT_MOD

