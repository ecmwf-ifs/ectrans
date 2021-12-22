      SUBROUTINE MXVA(SA,IAC,IAR,SB,IB,SC,IC,NRA,NCA)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: IAC
      INTEGER(KIND=JPIM) :: IAR
      INTEGER(KIND=JPIM) :: IB
      INTEGER(KIND=JPIM) :: IC
      INTEGER(KIND=JPIM) :: NRA
      INTEGER(KIND=JPIM) :: NCA
      INTEGER(KIND=JPIM) :: I,J
!RJ-warning F66 feature, overbounds
!RJ turning to assumed size until explicit interfaces (:)
      REAL(KIND=JPRB) :: SA(*)
      REAL(KIND=JPRB) :: SB(*)
      REAL(KIND=JPRB) :: SC(*)
!RJ       dimension sa(1),sb(1),sc(1)
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('MXVA',0,ZHOOK_HANDLE)
      IF(IC.EQ.1 .AND. IB.EQ.1) THEN
        SC(1:NRA)=0.0
        DO J=1,NCA
          DO I=1,NRA
            SC(I)=SC(I)+SA(1+(I-1)*IAC+(J-1)*IAR)*SB(J)
          ENDDO
        ENDDO
      ELSE
        DO I=1,NRA
          SC(1+(I-1)*IC)=0.0_JPRB
        ENDDO

        DO J=1,NCA
!OCL NOVREC
          DO I=1,NRA
            SC(1+(I-1)*IC)=SC(1+(I-1)*IC)                               &
     &       +SA(1+(I-1)*IAC+(J-1)*IAR)                                 &
     &       *SB(1+(J-1)*IB)
          ENDDO
        ENDDO
      ENDIF
      IF (LHOOK) CALL DR_HOOK('MXVA',1,ZHOOK_HANDLE)
      ENDSUBROUTINE MXVA
