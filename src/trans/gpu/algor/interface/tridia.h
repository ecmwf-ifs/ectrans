INTERFACE
SUBROUTINE TRIDIA(KN,KSYS,KFIRST,KEND,KTYP,PM,PRHS,PSOL)

! -----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSYS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIRST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PM(1+(KTYP-1)*(KSYS-1),KN,-1:1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHS(KSYS,KN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSOL(KSYS,KN)

! -----------------------------------------------------------------------------

END SUBROUTINE TRIDIA
END INTERFACE
