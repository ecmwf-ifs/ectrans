MODULE LTDIR_CTL_MOD
CONTAINS
SUBROUTINE LTDIR_CTL(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC)

!**** *LTDIR_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CTL(...)

!     Explicit arguments : 
!     -------------------- 
!     KF_FS      - number of fields in Fourier space
!     KF_UV      - local number of spectral u-v fields
!     KF_SCALARS - local number of scalar spectral fields
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
!     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR

USE LTDIR_MOD
USE TRLTOM_MOD

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)

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
IF(.NOT. ALLOCATED(FOUBUF))  ALLOCATE(FOUBUF(IBLEN))
CALL GSTATS(153,0)
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*KF_FS)
CALL GSTATS(153,1)
IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(FOUBUF_IN)
ENDIF

! Direct Legendre transform

ILED2 = 2*KF_FS
CALL GSTATS(103,0)
IF(KF_FS>0) THEN
#ifndef HLOMP
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
#endif
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTDIR(IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC)
  ENDDO
#ifndef HLOMP
!$OMP END PARALLEL DO
#endif
ENDIF
CALL GSTATS(103,1)

IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(FOUBUF)
ENDIF

!     -----------------------------------------------------------------

END SUBROUTINE LTDIR_CTL
END MODULE LTDIR_CTL_MOD
