      SUBROUTINE CTCAB_1DV (N,U,V)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK

      IMPLICIT NONE

      INTEGER(KIND=JPIM) :: N
      REAL(KIND=JPRB) :: U(N)
      REAL(KIND=JPRB) :: V(N)
!
      INTEGER(KIND=JPIM) :: I
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CTCAB_1DV',0,ZHOOK_HANDLE)
      DO 1 I=1,N
          V(I)=U(I)
 1    CONTINUE
      IF (LHOOK) CALL DR_HOOK('CTCAB_1DV',1,ZHOOK_HANDLE)
      ENDSUBROUTINE CTCAB_1DV
!
!-----------------------------------------------------------------------
!
