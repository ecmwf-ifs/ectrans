MODULE LTDIR_CTLAD_MOD
CONTAINS
SUBROUTINE LTDIR_CTLAD(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC)

!**** *LTDIR_CTLAD* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CTLAD(...)

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

USE LTDIRAD_MOD
USE TRMTOL_MOD

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
INTEGER_M,OPTIONAL,INTENT(IN)   :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN)   :: KFLDPTRSC(:)

INTEGER_M :: JM,IM,IBLEN,ILED2

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

IBLEN = D%NLENGT0B*2*KF_FS
IF(LALLOPERM) THEN
  IF(ALLOCATED(FOUBUF)) THEN
    IF(SIZE(FOUBUF) < IBLEN) THEN
      DEALLOCATE(FOUBUF)
    ENDIF
  ENDIF
ENDIF
IF(.NOT. ALLOCATED(FOUBUF)) ALLOCATE(FOUBUF(IBLEN))


! Direct Legendre transform

ILED2 = 2*KF_FS
CALL GSTATS(131,0)
IF(KF_FS > 0) THEN
#ifndef HLOMP
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JM,IM)
#endif
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTDIRAD(IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC)
  ENDDO
#ifndef HLOMP
!$OMP END PARALLEL DO
#endif
ENDIF
CALL GSTATS(131,1)

IBLEN = D%NLENGT0B*2*KF_FS
IF(LALLOPERM) THEN
  IF(ALLOCATED(FOUBUF_IN)) THEN
    IF(SIZE(FOUBUF_IN) < IBLEN) THEN
      DEALLOCATE(FOUBUF_IN)
    ENDIF
  ENDIF
ENDIF
IF(.NOT.ALLOCATED(FOUBUF_IN)) ALLOCATE(FOUBUF_IN(IBLEN))

CALL GSTATS(181,0)
CALL TRMTOL(FOUBUF,FOUBUF_IN,2*KF_FS)
CALL GSTATS(181,1)
IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(FOUBUF)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LTDIR_CTLAD
END MODULE LTDIR_CTLAD_MOD
