      SUBROUTINE MLIS0 (YDGEOMETRY,YDFIELDS,YDMTRAJ,N,SIMUL,PROSCA,     &
     &                  XN,FN,FPN,T,TMIN,TMAX,D,G,AMD,AMF,IMP,IO,       &
     &                  LOGIC,NAP,NAPMAX,X,NITER,ZTEMP,YDVARBC)
!AUTOPROMOTE
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
      USE GEOMETRY_MOD , ONLY : GEOMETRY
      USE FIELDS_MOD , ONLY : FIELDS
      USE MTRAJ_MOD  , ONLY : MTRAJ
      USE VARBC_CLASS,ONLY: CLASS_VARBC
      USE CONTROL_VECTORS_MOD
!
      IMPLICIT NONE
!
      TYPE(GEOMETRY),    INTENT(INOUT) :: YDGEOMETRY
      TYPE(FIELDS),      INTENT(INOUT) :: YDFIELDS
      TYPE(MTRAJ),       INTENT(INOUT) :: YDMTRAJ
      TYPE(CLASS_VARBC), INTENT(INOUT) :: YDVARBC
!
! --- arguments
!
      EXTERNAL :: SIMUL,PROSCA
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: IMP
      INTEGER(KIND=JPIM) :: IO
      INTEGER(KIND=JPIM) :: LOGIC
      INTEGER(KIND=JPIM) :: NAP
      INTEGER(KIND=JPIM) :: NAPMAX
      INTEGER(KIND=JPIM) :: NITER
      REAL(KIND=JPRB) :: FN
      REAL(KIND=JPRB) :: FPN
      REAL(KIND=JPRB) :: T
      REAL(KIND=JPRB) :: TMIN
      REAL(KIND=JPRB) :: TMAX
      REAL(KIND=JPRB) :: AMD
      REAL(KIND=JPRB) :: AMF
      TYPE (CONTROL_VECTOR) :: XN
      TYPE (CONTROL_VECTOR) :: D
      TYPE (CONTROL_VECTOR) :: G
      TYPE (CONTROL_VECTOR) :: X
!
! --- variables locales
!
      INTEGER(KIND=JPIM) :: I,INDIC,INDICA,INDICD
      REAL(KIND=JPRB) :: TESF,TESD,TG,FG,FPG,TD,TA
      REAL(KIND=JPRB) :: FA,FPA,D2,F,FP,FFN,FD,FPD
      REAL(KIND=JPRB) :: Z,TEST,BARMIN,BARMUL,BARMAX,BARR
      REAL(KIND=JPRB) :: GAUCHE,DROITE,TAA,PS,ZMAXP
      TYPE (CONTROL_VECTOR) :: ZTEMP

!

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('MLIS0',0,ZHOOK_HANDLE)
      IF (.NOT.(N.GT.0 .AND. FPN.LT.0.0_JPRB .AND. T.GT.0.0_JPRB        &
     & .AND. TMAX.GT.0.0_JPRB .AND. AMF.GT.0.0_JPRB                     &
     & .AND. AMD.GT.AMF .AND. AMD.LT.1.0_JPRB)) THEN
        LOGIC=6
        GOTO 999
      ENDIF

      TESF=AMF*FPN
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
      CALL PROSCA (N,D,D,PS)
      D2=PS
!
!               elimination d'un t initial ridiculement petit
!
      IF(T.LE.TMIN) THEN
        T=TMIN
        IF(T.GT.TMAX) THEN
          IF (IMP.GT.0) WRITE (IO,1007)
          TMIN=TMAX
        ENDIF
      ENDIF

      DO
        IF (FN+T*FPN.LT.FN+0.9_JPRB*T*FPN) EXIT
        T=2.0_JPRB*T
      ENDDO

      INDICA=1
      LOGIC=0
      IF (T.GT.TMAX) THEN
          T=TMAX
          LOGIC=1
      ENDIF
      IF (IMP.GE.4) WRITE (IO,1000) FPN,D2,TMIN,TMAX
!
!     --- nouveau x
!
      X%DATA = XN%DATA + T * D%DATA
!
! --- boucle
!
      BOUCLE:DO
      NAP=NAP+1
      IF(NAP.GT.NAPMAX) THEN
          LOGIC=4
          FN=FG
          XN%DATA = XN%DATA + TG * D%DATA
          EXIT
      ENDIF
      INDIC=4
!
!     --- appel simulateur
!
      CALL SIMUL(YDGEOMETRY,YDFIELDS,YDMTRAJ,INDIC,N,X,F,G,             &
     &           NITER,YDVARBC)
      IF(INDIC.EQ.0) THEN
!
!         --- arret demande par l'utilisateur
!
          LOGIC=5
          FN=F
          XN = X
         EXIT
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
      CALL PROSCA (N,D,G,PS)
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
  320 CONTINUE
      FN=F
      XN =X
      EXIT
!
!
!
  350 CONTINUE
      TG=T
      FG=F
      FPG=FP
      IF(TD.EQ.0.) THEN
!
!              extrapolation
!
        TAA=T
        GAUCHE=(1.0_JPRB+BARMIN)*T
        DROITE=10.0_JPRB*T
        CALL ECUBE (T,F,FP,TA,FA,FPA,GAUCHE,DROITE)
        TA=TAA
        IF(T.GE.TMAX) THEN
          LOGIC=1
          T=TMAX
        ENDIF
        GO TO 900
      ENDIF
!
!              interpolation
!
 500  CONTINUE
      IF(INDICA.LE.0) THEN
          TA=T
          T=0.9_JPRB*TG+0.1_JPRB*TD
      ELSE
        TEST=BARR*(TD-TG)
        GAUCHE=TG+TEST
        DROITE=TD-TEST
        TAA=T
        CALL ECUBE (T,F,FP,TA,FA,FPA,GAUCHE,DROITE)
        TA=TAA
        IF (T.GT.GAUCHE .AND. T.LT.DROITE) THEN
          BARR=MAX(BARMIN,BARR/BARMUL)
!         barr=barmin
        ELSE
          BARR=MIN(BARMUL*BARR,BARMAX)
        ENDIF
      ENDIF
!
! --- fin de boucle
!     - t peut etre bloque sur tmax
!       (venant de l'extrapolation avec logic=1)
!
  900 CONTINUE
      FA=F
      FPA=FP
  905 CONTINUE
      INDICA=INDIC
!
! --- faut-il continuer ?
!
      IF (TD.NE.0.0_JPRB) THEN
      IF (TD-TG.GE.TMIN) THEN
!
!     --- limite de precision machine (arret de secours) ?
!
      ZTEMP=0.0_JPRB
!OCL NOALIAS
      DO I=1,D%NSIZEL
          Z=XN%DATA(I)+T*D%DATA(I)
          ZTEMP%DATA(I)=MAX(ABS(Z-XN%DATA(I)),ABS(Z-X%DATA(I)))
      ENDDO
      ZMAXP=MAXVAL(ZTEMP)

      IF(ZMAXP.NE.0.) GO TO 950
!
! --- arret sur dxmin ou de secours
!
      ENDIF
      LOGIC=6
!
!     si indicd<0, derniers calculs non faits par simul
!
      IF (INDICD.LT.0) LOGIC=INDICD
!
!     si tg=0, xn = xn_depart,
!     sinon on prend xn=x_gauche qui fait decroitre f
!
      IF (TG.NE.0.0_JPRB) THEN
      FN=FG
      XN%DATA = XN%DATA + TG * D%DATA

      ENDIF
      IF (IMP.LE.0) EXIT
      WRITE (IO,1001)
      WRITE (IO,1005) TG,FG,FPG
      IF (LOGIC.EQ.6) WRITE (IO,1005) TD,FD,FPD
      IF (LOGIC.EQ.7) WRITE (IO,1006) TD,INDICD
      EXIT
!
!               recopiage de x et boucle
!
      ENDIF
  950 CONTINUE
      X%DATA = XN%DATA + T * D%DATA

      ENDDO BOUCLE

  999 CONTINUE
      IF (LHOOK) CALL DR_HOOK('MLIS0',1,ZHOOK_HANDLE)
      RETURN

 1000 FORMAT (/4X," mlis0   ",4X,"fpn=",E10.3," d2=",E9.2,              &
     & "  tmin=",E9.2," tmax=",E9.2)
 1001 FORMAT (/4X," mlis0",3X,"stop on tmin",8X,                        &
     &   "step",11X,"functions",5X,"derivatives")
 1002 FORMAT (4X," mlis0",37X,E10.3,2E11.3)
 1003 FORMAT (4X," mlis0",E14.3,2E11.3)
 1004 FORMAT (4X," mlis0",37X,E10.3," indic=",I3)
 1005 FORMAT (4X," mlis0",14X,2E18.8,E11.3)
 1006 FORMAT (4X," mlis0",14X,E18.8,"      indic=",I3)
 1007 FORMAT (/4X," mlis0",10X,"tmin forced to tmax")
 1008 FORMAT (/4X," mlis0",10X,"inconsistent call")

      ENDSUBROUTINE MLIS0
!
!-----------------------------------------------------------------------
!
