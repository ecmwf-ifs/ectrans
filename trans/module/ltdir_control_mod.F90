MODULE LTDIR_CONTROL_MOD
CONTAINS
SUBROUTINE LTDIR_CONTROL(PSPVOR,PSPDIV,PSPSCALAR)

!**** *LTDIR_CONTROL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CONTROL(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR

USE LTDIR_MOD
USE TRLTOM_MOD

IMPLICIT NONE

REAL_B ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)

INTEGER_M :: JM,IM

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

ALLOCATE(FOUBUF(D%NLENGT0B*2*NF_FS))
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*NF_FS)
DEALLOCATE(FOUBUF_IN)

! Direct Legendre transform

NLED2 = 2*NF_FS
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JM,IM)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)
  CALL LTDIR(IM,JM,PSPVOR,PSPDIV,PSPSCALAR)
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(FOUBUF)

!     ------------------------------------------------------------------

END SUBROUTINE LTDIR_CONTROL
END MODULE LTDIR_CONTROL_MOD
