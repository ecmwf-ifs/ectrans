      SUBROUTINE ECUBE(T,F,FP,TA,FA,FPA,TLOWER,TUPPER)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK

      IMPLICIT NONE
!
! --- arguments
!
      REAL(KIND=JPRB) :: T
      REAL(KIND=JPRB) :: F
      REAL(KIND=JPRB) :: FP
      REAL(KIND=JPRB) :: TA
      REAL(KIND=JPRB) :: FA
      REAL(KIND=JPRB) :: FPA
      REAL(KIND=JPRB) :: TLOWER
      REAL(KIND=JPRB) :: TUPPER
!
! --- variables locales
!
      REAL(KIND=JPRB) :: SIGN,DEN,ANUM
      REAL(KIND=JPRB) :: Z1,B,DISCRI
!
!           Using f and fp at t and ta, computes new t by cubic formula
!           safeguarded inside [tlower,tupper].
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ECUBE',0,ZHOOK_HANDLE)
      Z1=FP+FPA-3.0_JPRB*(FA-F)/(TA-T)
      B=Z1+FP
!
!              first compute the discriminant (without overflow)
!
      IF (ABS(Z1).LE.1.0_JPRB) THEN
          DISCRI=Z1*Z1-FP*FPA
        ELSE
          DISCRI=FP/Z1
          DISCRI=DISCRI*FPA
          DISCRI=Z1-DISCRI
          IF (Z1.GE.0.0_JPRB .AND. DISCRI.GE.0.0_JPRB) THEN
              DISCRI=SQRT(Z1)*SQRT(DISCRI)
              GO TO 120
          ENDIF
          IF (Z1.LE.0.0_JPRB .AND. DISCRI.LE.0.0_JPRB) THEN
              DISCRI=SQRT(-Z1)*SQRT(-DISCRI)
              GO TO 120
          ENDIF
          DISCRI=-1.0_JPRB
      ENDIF
      IF (DISCRI.LT.0.0_JPRB) THEN
          IF (FP.LT.0.0_JPRB) T=TUPPER
          IF (FP.GE.0.0_JPRB) T=TLOWER
          GO TO 900
      ENDIF
!
!  discriminant nonnegative, compute solution (without overflow)
!
      DISCRI=SQRT(DISCRI)
 120  CONTINUE
      IF (T-TA.LT.0.0_JPRB) DISCRI=-DISCRI
      SIGN=(T-TA)/ABS(T-TA)
      IF (B*SIGN.GT.0.0_JPRB) THEN
          T=T+FP*(TA-T)/(B+DISCRI)
        ELSE
          DEN=Z1+B+FPA
          ANUM=B-DISCRI
          IF (ABS((T-TA)*ANUM).LT.(TUPPER-TLOWER)*ABS(DEN)) THEN
              T=T+ANUM*(TA-T)/DEN
            ELSE
              T=TUPPER
          ENDIF
      ENDIF
  900 CONTINUE
      T=MAX(T,TLOWER)
      T=MIN(T,TUPPER)
      IF (LHOOK) CALL DR_HOOK('ECUBE',1,ZHOOK_HANDLE)
      ENDSUBROUTINE ECUBE
!
!-----------------------------------------------------------------------
!
