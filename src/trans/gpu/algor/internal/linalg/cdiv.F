      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
      IMPLICIT NONE
!
      REAL(KIND=JPRB) :: AR
      REAL(KIND=JPRB) :: AI
      REAL(KIND=JPRB) :: BR
      REAL(KIND=JPRB) :: BI
      REAL(KIND=JPRB) :: CR
      REAL(KIND=JPRB) :: CI
!
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
!
      REAL(KIND=JPRB) :: S,ARS,AIS,BRS,BIS
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CDIV',0,ZHOOK_HANDLE)
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      IF (LHOOK) CALL DR_HOOK('CDIV',1,ZHOOK_HANDLE)
      ENDSUBROUTINE CDIV
