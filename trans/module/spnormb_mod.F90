MODULE SPNORMB_MOD
CONTAINS
SUBROUTINE SPNORMB(PSPEC,KFLD,PMET,PSM)

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_DISTR

REAL_B    ,INTENT(IN)  :: PSPEC(:,:)
REAL_B    ,INTENT(IN)  :: PMET(0:R%NSMAX)
INTEGER_M ,INTENT(IN)  :: KFLD
REAL_B    ,INTENT(OUT) :: PSM(:,:)

INTEGER_M :: JM ,JFLD ,JN ,IM ,ISP

!     ------------------------------------------------------------------

PSM(:,:) = _ZERO_

!$OMP PARALLEL DO SCHEDULE(STATIC,1)  PRIVATE(JN,JFLD)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)
  IF(IM == 0)THEN
    DO JN=0,R%NSMAX
      ISP = D%NASM0(0)+JN*2
      DO JFLD=1,KFLD
        PSM(JFLD,JM) = PSM(JFLD,JM)+PMET(JN)*PSPEC(JFLD,ISP)**2
      ENDDO
    ENDDO
  ELSE
    DO JN=IM,R%NSMAX
      ISP = D%NASM0(IM)+(JN-IM)*2
      DO JFLD=1,KFLD
        PSM(JFLD,JM) = PSM(JFLD,JM)+_TWO_*PMET(JN)*&
         &(PSPEC(JFLD,ISP)**2+PSPEC(JFLD,ISP+1)**2)
      ENDDO
    ENDDO
  ENDIF
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

END SUBROUTINE SPNORMB
END MODULE SPNORMB_MOD





