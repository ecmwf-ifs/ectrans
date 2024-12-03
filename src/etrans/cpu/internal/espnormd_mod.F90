! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ESPNORMD_MOD
CONTAINS
SUBROUTINE ESPNORMD(PSPEC,KFLD,PMET,PSM)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D

USE TPMALD_DISTR    ,ONLY : DALD
!

IMPLICIT NONE

REAL(KIND=JPRB)    ,INTENT(IN)  :: PSPEC(:,:)
REAL(KIND=JPRB)    ,INTENT(IN)  :: PMET(0:R%NSPEC_G)
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KFLD
REAL(KIND=JPRB)    ,INTENT(OUT) :: PSM(:,:)
INTEGER(KIND=JPIM) :: JM ,JFLD ,JN ,IM ,ISP
INTEGER(KIND=JPIM) :: IN,ISPE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESPNORMD_MOD:ESPNORMD',0,ZHOOK_HANDLE)

!$OMP PARALLEL DO SCHEDULE(STATIC,1)  PRIVATE(JM,IM,JN,ISP,JFLD,IN,ISPE)
DO JM=1,D%NUMP
  PSM(:,JM) = 0.0_JPRB
  IM = D%MYMS(JM)

  IN=DALD%NCPL2M(IM)/2 - 1
  DO JN=0,IN
    ISP=DALD%NESM0(IM) + (JN)*4
    ISPE=DALD%NPME (IM) + JN
    DO JFLD=1,KFLD
      PSM(JFLD,JM) =PSM(JFLD,JM)&
       & + PMET(ISPE) *&
       & ( PSPEC(JFLD,ISP  )**2 + PSPEC(JFLD,ISP+1)**2 +&
       & PSPEC(JFLD,ISP+2)**2 + PSPEC(JFLD,ISP+3)**2   )

    ENDDO
  ENDDO

ENDDO
!$OMP END PARALLEL DO

IF (LHOOK) CALL DR_HOOK('ESPNORMD_MOD:ESPNORMD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ESPNORMD
END MODULE ESPNORMD_MOD

