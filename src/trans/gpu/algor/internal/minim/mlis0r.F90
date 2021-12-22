      SUBROUTINE MLIS0R (YDGEOMETRY,                                    &
     &                   N,SIMULR,PROSCAR,XN,FN,FPN,T,TMIN,TMAX,D,G,    &
     &                   AMD,AMF,IMP,IO,LOGIC,NAP,NAPMAX,X,NITER)
!AUTOPROMOTE

      USE GEOMETRY_MOD , ONLY : GEOMETRY
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
! ----
!
!     mlis0 + minuscules + commentaires
!     + version amelioree (XII 88): interpolation cubique systematique
!       et anti-overflows
!     + declaration variables (II/89, JCG).
!     + barr is also progressively decreased (12/93, CL & JChG).
!       barmul is set to 5.
!
!     ----------------------------------------------------------------
!
!        en sortie logic =
!
!        0          descente serieuse
!        1          descente bloquee
!        4          nap > napmax
!        5          retour a l'utilisateur
!        6          fonction et gradient pas d'accord
!        < 0        contrainte implicite active
!
! ----
!
! --- arguments
!
      IMPLICIT NONE
!
      TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
      EXTERNAL :: SIMULR,PROSCAR
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: IMP
      INTEGER(KIND=JPIM) :: IO
      INTEGER(KIND=JPIM) :: LOGIC
      INTEGER(KIND=JPIM) :: NAP
      INTEGER(KIND=JPIM) :: NAPMAX
      INTEGER(KIND=JPIM) :: NITER
      REAL(KIND=JPRB) :: XN(N)
      REAL(KIND=JPRB) :: FN
      REAL(KIND=JPRB) :: FPN
      REAL(KIND=JPRB) :: T
      REAL(KIND=JPRB) :: TMIN
      REAL(KIND=JPRB) :: TMAX
      REAL(KIND=JPRB) :: D(N)
      REAL(KIND=JPRB) :: G(N)
      REAL(KIND=JPRB) :: AMD
      REAL(KIND=JPRB) :: AMF
      REAL(KIND=JPRB) :: X(N)
!
! --- variables locales
!
      INTEGER(KIND=JPIM) :: I,INDIC,INDICA,INDICD
      REAL(KIND=JPRB) :: TESF,TESD,TG,FG,FPG,TD,TA
      REAL(KIND=JPRB) :: FA,FPA,D2,F,FP,FFN,FD,FPD
      REAL(KIND=JPRB) :: Z,TEST,BARMIN,BARMUL,BARMAX,BARR
      REAL(KIND=JPRB) :: GAUCHE,DROITE,TAA
      REAL(KIND=JPRB) :: PS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
 1000 FORMAT (/4X,9H MLIS0   ,4X,4HFPN=,E10.3,4H D2=,E9.2,              &
     & 7H  TMIN=,E9.2,6H TMAX=,E9.2)
 1001 FORMAT (/4X,7H MLIS0R,3X,"stop on tmin",8X,                       &
     &   "step",11X,"functions",5X,"derivatives")
 1002 FORMAT (4X,7H MLIS0R,37X,E10.3,2E11.3)
 1003 FORMAT (4X,7H MLIS0R,E14.3,2E11.3)
 1004 FORMAT (4X,7H MLIS0R,37X,E10.3,7H INDIC=,I3)
 1005 FORMAT (4X,7H MLIS0R,14X,2E18.8,E11.3)
 1006 FORMAT (4X,7H MLIS0R,14X,E18.8,12H      INDIC=,I3)
 1007 FORMAT (/4X,7H MLIS0R,10X,"tmin forced to tmax")
 1008 FORMAT (/4X,7H MLIS0R,10X,"inconsistent call")
      IF (LHOOK) CALL DR_HOOK('MLIS0R',0,ZHOOK_HANDLE)
      ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM)
      IF (N.GT.0 .AND. FPN.LT.0. .AND. T.GT.0.                          &
     & .AND. TMAX.GT.0.0_JPRB .AND. AMF.GT.0.0_JPRB                     &
     & .AND. AMD.GT.AMF .AND. AMD.LT.1.0_JPRB) GO TO 5
      LOGIC=6
      GO TO 999
    5 TESF=AMF*FPN
      TESD=AMD*FPN
      BARMIN=0.01_JPRB
      BARMUL=5.0_JPRB
      BARMAX=0.3_JPRB
      BARR=BARMIN
      TD=0.0_JPRB
      TG=0.0_JPRB
      FG=FN
      FPG=FPN
      TA=0.0_JPRB
      FA=FN
      FPA=FPN
      CALL PROSCAR (YDDIM,N,D,D,PS)
      D2=PS
!
!               elimination d'un t initial ridiculement petit
!
      IF (T.GT.TMIN) GO TO 20
      T=TMIN
      IF (T.LE.TMAX) GO TO 20
      IF (IMP.GT.0) WRITE (IO,1007)
      TMIN=TMAX
   20 IF (FN+T*FPN.LT.FN+0.9_JPRB*T*FPN) GO TO 30
      T=2.*T
      GO TO 20
   30 INDICA=1
      LOGIC=0
      IF (T.GT.TMAX) THEN
          T=TMAX
          LOGIC=1
      ENDIF
      IF (IMP.GE.4) WRITE (IO,1000) FPN,D2,TMIN,TMAX
!
!     --- nouveau x
!
      DO 50 I=1,N
          X(I)=XN(I)+T*D(I)
   50 CONTINUE
!
! --- boucle
!
  100 NAP=NAP+1
      IF(NAP.GT.NAPMAX) THEN
          LOGIC=4
          FN=FG
          DO 120 I=1,N
              XN(I)=XN(I)+TG*D(I)
  120     CONTINUE
          GO TO 999
      ENDIF
      INDIC=4
!
!     --- appel simulateur
!
      CALL SIMULR(YDGEOMETRY,INDIC,N,X,F,G,NITER)
      IF(INDIC.EQ.0) THEN
!
!         --- arret demande par l'utilisateur
!
          LOGIC=5
          FN=F
          DO 170 I=1,N
              XN(I)=X(I)
  170     CONTINUE
          GO TO 999
      ENDIF
      IF(INDIC.LT.0) THEN
!
!         --- les calculs n'ont pas pu etre effectues par le simulateur
!
          TD=T
          INDICD=INDIC
          LOGIC=0
          IF (IMP.GE.4) WRITE (IO,1004) T,INDIC
          T=TG+0.1_JPRB*(TD-TG)
          GO TO 905
      ENDIF
!
!     --- les tests elementaires sont faits, on y va
!
      CALL PROSCAR (YDDIM,N,D,G,PS)
      FP=PS
!
!     --- premier test de Wolfe
!
      FFN=F-FN
      IF(FFN.GT.T*TESF) THEN
          TD=T
          FD=F
          FPD=FP
          INDICD=INDIC
          LOGIC=0
          IF(IMP.GE.4) WRITE (IO,1002) T,FFN,FP
          GO TO 500
      ENDIF
!
!     --- test 1 ok, donc deuxieme test de Wolfe
!
      IF(IMP.GE.4) WRITE (IO,1003) T,FFN,FP
      IF(FP.GT.TESD) THEN
          LOGIC=0
          GO TO 320
      ENDIF
      IF (LOGIC.EQ.0) GO TO 350
!
!     --- test 2 ok, donc pas serieux, on sort
!
  320 FN=F
      DO 330 I=1,N
          XN(I)=X(I)
  330 CONTINUE
      GO TO 999
!
!
!
  350 TG=T
      FG=F
      FPG=FP
      IF(TD.NE.0.0_JPRB) GO TO 500
!
!              extrapolation
!
      TAA=T
      GAUCHE=(1.0_JPRB+BARMIN)*T
      DROITE=10.0_JPRB*T
      CALL ECUBER (T,F,FP,TA,FA,FPA,GAUCHE,DROITE)
      TA=TAA
      IF(T.LT.TMAX) GO TO 900
      LOGIC=1
      T=TMAX
      GO TO 900
!
!              interpolation
!
  500 IF(INDICA.LE.0) THEN
          TA=T
          T=0.9_JPRB*TG+0.1_JPRB*TD
          GO TO 900
      ENDIF
      TEST=BARR*(TD-TG)
      GAUCHE=TG+TEST
      DROITE=TD-TEST
      TAA=T
      CALL ECUBER (T,F,FP,TA,FA,FPA,GAUCHE,DROITE)
      TA=TAA
      IF (T.GT.GAUCHE .AND. T.LT.DROITE) THEN
          BARR=MAX(BARMIN,BARR/BARMUL)
!         barr=barmin
        ELSE
          BARR=MIN(BARMUL*BARR,BARMAX)
      ENDIF
!
! --- fin de boucle
!     - t peut etre bloque sur tmax
!       (venant de l'extrapolation avec logic=1)
!
  900 FA=F
      FPA=FP
  905 INDICA=INDIC
!
! --- faut-il continuer ?
!
      IF (TD.EQ.0.0_JPRB) GO TO 950
      IF (TD-TG.LT.TMIN) GO TO 920
!
!     --- limite de precision machine (arret de secours) ?
!
      DO 910 I=1,N
          Z=XN(I)+T*D(I)
          IF (Z.NE.XN(I).AND.Z.NE.X(I)) GO TO 950
  910 CONTINUE
!
! --- arret sur dxmin ou de secours
!
  920 LOGIC=6
!
!     si indicd<0, derniers calculs non faits par simul
!
      IF (INDICD.LT.0) LOGIC=INDICD
!
!     si tg=0, xn = xn_depart,
!     sinon on prend xn=x_gauche qui fait decroitre f
!
      IF (TG.EQ.0.0_JPRB) GO TO 940
      FN=FG
      DO 930 I=1,N
  930 XN(I)=XN(I)+TG*D(I)
  940 IF (IMP.LE.0) GO TO 999
      WRITE (IO,1001)
      WRITE (IO,1005) TG,FG,FPG
      IF (LOGIC.EQ.6) WRITE (IO,1005) TD,FD,FPD
      IF (LOGIC.EQ.7) WRITE (IO,1006) TD,INDICD
      GO TO 999
!
!               recopiage de x et boucle
!
  950 DO 960 I=1,N
  960 X(I)=XN(I)+T*D(I)
      GO TO 100
  999 CONTINUE
      END ASSOCIATE
      IF (LHOOK) CALL DR_HOOK('MLIS0R',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE MLIS0R
!
!-----------------------------------------------------------------------
!
