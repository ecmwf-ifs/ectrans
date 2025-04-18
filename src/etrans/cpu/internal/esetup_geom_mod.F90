! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ESETUP_GEOM_MOD
CONTAINS
SUBROUTINE ESETUP_GEOM

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV
USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE TPMALD_DIM      ,ONLY : RALD
!USE TPM_FIELDS
USE TPM_GEOMETRY    ,ONLY : G
!

IMPLICIT NONE

INTEGER(KIND=JPIM) :: IDGLU(0:RALD%NMSMAX,R%NDGNH)
INTEGER(KIND=JPIM) :: JGL,JM

LOGICAL    :: LLP1,LLP2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESETUP_GEOM_MOD:ESETUP_GEOM',0,ZHOOK_HANDLE)
IF(.NOT.D%LGRIDONLY) THEN
LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1

IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_GEOM ==='

ALLOCATE (G%NMEN(R%NDGL))
IF(LLP2)WRITE(NOUT,9) 'G%NMEN   ',SIZE(G%NMEN   ),SHAPE(G%NMEN   )
G%NMEN(:)=RALD%NMSMAX
IF(LLP1) THEN
  WRITE(NOUT,FMT='('' (JGL,G%NLOEN,G%NMEN) '')')
  WRITE(NOUT,FMT='(8(1X,''('',I4,I4,I4,'')''))')&
   & (JGL,G%NLOEN(JGL),G%NMEN(JGL),JGL=1,R%NDGL)
ENDIF
ALLOCATE(G%NDGLU(0:RALD%NMSMAX))
IF(LLP2)WRITE(NOUT,9) 'G%NDGLU   ',SIZE(G%NDGLU   ),SHAPE(G%NDGLU   )
IDGLU(:,:) = 0
G%NDGLU(:) = 0
DO JGL=1,R%NDGNH
  DO JM=0,G%NMEN(JGL)
    IDGLU(JM,JGL) = 1
  ENDDO
ENDDO
DO JM=0,RALD%NMSMAX
  DO JGL=1,R%NDGNH
    G%NDGLU(JM) = G%NDGLU(JM)+IDGLU(JM,JGL)
  ENDDO
ENDDO
IF(LLP1) THEN
  WRITE(NOUT,FMT='('' (JM,G%NDGLU) '')')
  WRITE(NOUT,FMT='(10(1X,''('',I4,I4,'')''))')&
   & (JM,G%NDGLU(JM),JM=0,RALD%NMSMAX)
ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('ESETUP_GEOM_MOD:ESETUP_GEOM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE ESETUP_GEOM
END MODULE ESETUP_GEOM_MOD
