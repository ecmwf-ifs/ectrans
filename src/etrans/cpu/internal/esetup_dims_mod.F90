! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ESETUP_DIMS_MOD
CONTAINS
SUBROUTINE ESETUP_DIMS

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPMALD_DIM      ,ONLY : RALD
!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: JM,JN,ISPOLEG
INTEGER(KIND=JPIM) :: ISMAX(0:R%NSMAX),ISNAX(0:RALD%NMSMAX)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESETUP_DIMS_MOD:ESETUP_DIMS',0,ZHOOK_HANDLE)
ISPOLEG = 0
DO JM=0,R%NSMAX
  DO JN=JM,R%NTMAX+1
    ISPOLEG = ISPOLEG+1
  ENDDO
ENDDO
R%NSPOLEG = ISPOLEG
CALL ELLIPS(R%NSMAX,RALD%NMSMAX,ISNAX,ISMAX)
R%NSPEC_G=0
DO JM=0,RALD%NMSMAX
  R%NSPEC_G=R%NSPEC_G+2*(ISNAX(JM)+1)
ENDDO
R%NSPEC2_G = R%NSPEC_G*2

R%NDGNH = (R%NDGL+1)/2

R%NLEI1 = R%NSMAX+4+MOD(R%NSMAX+4+1,2)
R%NLEI3 = R%NDGNH+MOD(R%NDGNH+2,2)

R%NLED3 = R%NTMAX+2+MOD(R%NTMAX+3,2)
R%NLED4 = R%NTMAX+3+MOD(R%NTMAX+4,2)
IF (LHOOK) CALL DR_HOOK('ESETUP_DIMS_MOD:ESETUP_DIMS',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ESETUP_DIMS
END MODULE ESETUP_DIMS_MOD
