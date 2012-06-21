MODULE LTDIR_CTLAD_MOD
CONTAINS
SUBROUTINE LTDIR_CTLAD(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR, &
 & PSPSC3A,PSPSC3B,PSPSC2, &
 & KFLDPTRUV,KFLDPTRSC)

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
USE YOMGSTATS, ONLY : LSYNCSTATS

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPSC3B(:,:,:)
REAL_B ,OPTIONAL, INTENT(INOUT) :: PSPSC2(:,:)
INTEGER_M,OPTIONAL,INTENT(IN)   :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN)   :: KFLDPTRSC(:)

INTEGER_M :: JM,IM,IBLEN,ILED2

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

IBLEN = D%NLENGT0B*2*KF_FS
ALLOCATE(FOUBUF_IN(IBLEN))
ALLOCATE(FOUBUF(IBLEN))
! Direct Legendre transform

ILED2 = 2*KF_FS
IF (.NOT.LSYNCSTATS) CALL GSTATS(1646,0)
IF(KF_FS > 0) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTDIRAD(IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC)
  ENDDO
!$OMP END PARALLEL DO
ENDIF
IF (.NOT.LSYNCSTATS) CALL GSTATS(1646,1)

CALL GSTATS(181,0)
CALL TRMTOL(FOUBUF,FOUBUF_IN,2*KF_FS)
CALL GSTATS(181,1)
DEALLOCATE(FOUBUF)

!     ------------------------------------------------------------------

END SUBROUTINE LTDIR_CTLAD
END MODULE LTDIR_CTLAD_MOD
