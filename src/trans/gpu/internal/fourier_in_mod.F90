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

MODULE FOURIER_IN_MOD
CONTAINS
SUBROUTINE FOURIER_IN(FOUBUF,PREEL_COMPLEX,KFIELDS)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

!     Explicit arguments :  PREEL_COMPLEX - local fourier/GP array
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

USE PARKIND_ECTRANS ,ONLY : JPIM, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, D_NSTAGTF,D_NSTAGT0B,D_NPNTGTB0,D_NPROCM,D_NPTRLS
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN, G_NLOEN
USE TPM_GEN         ,ONLY : NOUT
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
REAL(KIND=JPRBT), ALLOCATABLE, INTENT(INOUT) :: FOUBUF(:)
REAL(KIND=JPRBT), INTENT(OUT) :: PREEL_COMPLEX(:,:)
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

OFFSET_VAR=D_NPTRLS(MYSETW)

#ifdef ACCGPU
!$ACC DATA PRESENT(G_NLOEN,G_NMEN,D_NPROCM,D_NSTAGT0B,D_NPNTGTB0,FOUBUF,PREEL_COMPLEX,D_NSTAGTF) ASYNC(1)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(PRESENT,ALLOC:G_NMEN,D_NPROCM,D_NSTAGT0B,D_MSTABF,D_NPNTGTB0,FOUBUF,PREEL_COMPLEX,D_NSTAGTF)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(IGLG,IPROC,ISTA) DEFAULT(NONE) &
!$OMP&       SHARED(IBEG,IEND,IINC,KFIELDS,D_NPTRLS,MYSETW,G_NMEN, &
!$OMP&              D_NPROCM,D_NSTAGT0B,D_MSTABF,D_NPNTGTB0,PREEL_COMPLEX,D_NSTAGTF,FOUBUF)
#endif

#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(IGLG,IOFF_COMPLEX,JM,IPROC,ISTA) &
!$ACC&       FIRSTPRIVATE(IBEG,IEND,IINC,KFIELDS,OFFSET_VAR) ASYNC(1)
#endif
DO KGL=IBEG,IEND,IINC
  DO JF=1,KFIELDS
    IGLG = OFFSET_VAR+KGL-1
    IOFF_COMPLEX = (D_NSTAGTF(KGL)/2)+1

    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      IPROC = D_NPROCM(JM)
      ISTA  = (D_NSTAGT0B(IPROC)+D_NPNTGTB0(JM,KGL))*KFIELDS*2

      PREEL_COMPLEX(2*JF-1,JM+IOFF_COMPLEX) = FOUBUF(ISTA+2*JF-1)
      PREEL_COMPLEX(2*JF,  JM+IOFF_COMPLEX) = FOUBUF(ISTA+2*JF  )
    ENDDO
    !$ACC LOOP SEQ
    DO JM=G_NMEN(IGLG)+1,(G_NLOEN(IGLG)+4)/2-1
      ! Truncation (not sure what is the exact upper bound here...)
      ! Same is also in FSC for the new fields. I *think* it should be N/2+1 elements in total
      ! TODO: Make sure this is correct
      PREEL_COMPLEX(2*JF-1,JM+IOFF_COMPLEX) = 0._JPRBT
      PREEL_COMPLEX(2*JF,  JM+IOFF_COMPLEX) = 0._JPRBT
    ENDDO
  ENDDO
ENDDO
#ifdef ACCGPU
!$ACC END DATA
!$ACC WAIT(1)
#endif
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif

#ifdef OMPGPU
#endif
#ifdef ACCGPU
!$ACC EXIT DATA DELETE(FOUBUF)
#endif
DEALLOCATE(FOUBUF)

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_IN
END MODULE FOURIER_IN_MOD
