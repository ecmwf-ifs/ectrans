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
SUBROUTINE FOURIER_IN(FOUBUF,PREEL_COMPLEX,KF_CURRENT,KF_TOTAL)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

!     Explicit arguments :  PREEL_COMPLEX - local fourier/GP array
!     --------------------  KF_CURRENT - number of fields that are read (from Legendre space)
!                           KF_TOTAL - total fields in PREEL ("stride")
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

REAL(KIND=JPRBT), ALLOCATABLE, INTENT(INOUT) :: FOUBUF(:)
REAL(KIND=JPRBT), INTENT(OUT) :: PREEL_COMPLEX(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_CURRENT, KF_TOTAL

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,ISTA,OFFSET_VAR,IOFF_COMPLEX,KGL
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

!$ACC DATA PRESENT(G_NLOEN,G_NMEN,D_NPROCM,D_NSTAGT0B,D_NPNTGTB0,FOUBUF,PREEL_COMPLEX,D_NSTAGTF) ASYNC(1)

OFFSET_VAR=D_NPTRLS(MYSETW)
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IGLG,IOFF_COMPLEX,JM,IPROC,ISTA) DEFAULT(NONE) ASYNC(1)
DO KGL=IBEG,IEND,IINC
  DO JF=1,KF_CURRENT
    IGLG = OFFSET_VAR+KGL-1
    IOFF_COMPLEX = D_NSTAGTF(KGL)/2

    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      IPROC = D_NPROCM(JM)
      ISTA  = (D_NSTAGT0B(IPROC)+D_NPNTGTB0(JM,KGL))*KF_CURRENT*2

      PREEL_COMPLEX(2*JF-1+2*KF_TOTAL*(JM+IOFF_COMPLEX)) = FOUBUF(ISTA+2*JF-1)
      PREEL_COMPLEX(2*JF  +2*KF_TOTAL*(JM+IOFF_COMPLEX)) = FOUBUF(ISTA+2*JF  )
    ENDDO
    !$ACC LOOP SEQ
    DO JM=G_NMEN(IGLG)+1,(G_NLOEN(IGLG)+4)/2-1
      ! Truncation (not sure what is the exact upper bound here...)
      ! Same is also in FSC for the new fields. I *think* it should be N/2+1 elements in total
      ! TODO: Make sure this is correct
      PREEL_COMPLEX(2*JF-1+2*KF_TOTAL*(JM+IOFF_COMPLEX)) = 0._JPRBT
      PREEL_COMPLEX(2*JF  +2*KF_TOTAL*(JM+IOFF_COMPLEX)) = 0._JPRBT
    ENDDO
  ENDDO
ENDDO
!$ACC END DATA

!$ACC WAIT(1)

!$ACC EXIT DATA DELETE(FOUBUF)
DEALLOCATE(FOUBUF)

END SUBROUTINE FOURIER_IN
END MODULE FOURIER_IN_MOD

