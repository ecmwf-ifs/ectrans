MODULE LTDIR_CTL_MOD
CONTAINS
SUBROUTINE LTDIR_CTL(PSPVOR,PSPDIV,PSPSCALAR)

!**** *LTDIR_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CTL(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR

USE LTDIR_MOD
USE TRLTOM_MOD

IMPLICIT NONE

REAL_B ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)

INTEGER_M :: JM,IM,IBLEN

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

IBLEN = D%NLENGT0B*2*NF_FS
IF(LALLOPERM) THEN
  IF(ALLOCATED(FOUBUF)) THEN
    IF(SIZE(FOUBUF) < IBLEN) THEN
      DEALLOCATE(FOUBUF)
    ENDIF
  ENDIF
ENDIF
IF(.NOT. ALLOCATED(FOUBUF))  ALLOCATE(FOUBUF(IBLEN))
CALL GSTATS(153,0)
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*NF_FS)
CALL GSTATS(153,1)
IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(FOUBUF_IN)
ENDIF

! Direct Legendre transform

NLED2 = 2*NF_FS
CALL GSTATS(103,0)
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JM,IM)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)
  CALL LTDIR(IM,JM,PSPVOR,PSPDIV,PSPSCALAR)
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(103,1)

IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(FOUBUF)
ENDIF

!     -----------------------------------------------------------------

END SUBROUTINE LTDIR_CTL
END MODULE LTDIR_CTL_MOD
