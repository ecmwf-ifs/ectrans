      SUBROUTINE DDR (YDDIM,PROSCAR,CTONBR,CTCABR,N,SSCALE,NM,DEPL,AUX, &
     &                JMIN,JMAX,PRECOS,DIAG,ALPHA,YBAR,SBAR)
!AUTOPROMOTE

      USE YOMDIM  , ONLY : TDIM
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!----
!
!     calcule le produit H.g ou
!         . H est une matrice construite par la formule de bfgs inverse
!           a nm memoires a partir de la matrice diagonale diag
!           dans un espace hilbertien dont le produit scalaire
!           est donne par proscar
!           (cf. J. Nocedal, Math. of Comp. 35/151 (1980) 773-782)
!         . g est un vecteur de dimension n (en general le gradient)
!
!     la matrice diag apparait donc comme un preconditionneur diagonal
!
!     depl = g (en entree), = H g (en sortie)
!
!     la matrice H est memorisee par les vecteurs des tableaux
!     ybar, sbar et les pointeurs jmin, jmax
!
!     alpha(nm) est une zone de travail
!
!     izs(1),rzs(1),dzs(1) sont des zones de travail pour proscar
!
!----
!
!         arguments
!
      TYPE(TDIM) , INTENT(INOUT) :: YDDIM
      LOGICAL :: SSCALE
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: NM
      INTEGER(KIND=JPIM) :: JMIN
      INTEGER(KIND=JPIM) :: JMAX
      REAL(KIND=JPRB) :: DEPL(N)
      REAL(KIND=JPRB) :: PRECOS
      REAL(KIND=JPRB) :: DIAG(N)
      REAL(KIND=JPRB) :: ALPHA(NM)
      REAL(KIND=JPRB) :: YBAR(N,1)
      REAL(KIND=JPRB) :: SBAR(N,1)
      REAL(KIND=JPRB) :: AUX(N)
      EXTERNAL :: PROSCAR,CTONBR,CTCABR
!
!         variables locales
!
      INTEGER(KIND=JPIM) :: JFIN,I,J,JP
      REAL(KIND=JPRB) :: R
      REAL(KIND=JPRB) :: PS
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('DDR',0,ZHOOK_HANDLE)
      JFIN=JMAX
      IF (JFIN.LT.JMIN) JFIN=JMAX+NM
!
!         phase de descente
!
      DO 100 J=JFIN,JMIN,-1
          JP=J
          IF (JP.GT.NM) JP=JP-NM
          CALL PROSCAR (YDDIM,N,DEPL,SBAR(1,JP),PS)
          R=PS
          ALPHA(JP)=R
          DO 20 I=1,N
              DEPL(I)=DEPL(I)-R*YBAR(I,JP)
20        CONTINUE
100   CONTINUE
!
!         preconditionnement
!
      IF (SSCALE) THEN
          DO 150 I=1,N
              DEPL(I)=DEPL(I)*PRECOS
  150     CONTINUE
      ELSE
          CALL CTONBR (N,DEPL,AUX)
          DO 151 I=1,N
              AUX(I)=AUX(I)*DIAG(I)
  151     CONTINUE
          CALL CTCABR (N,AUX,DEPL)
      ENDIF
!
!         remontee
!
      DO 200 J=JMIN,JFIN
          JP=J
          IF (JP.GT.NM) JP=JP-NM
          CALL PROSCAR (YDDIM,N,DEPL,YBAR(1,JP),PS)
          R=ALPHA(JP)-PS
          DO 120 I=1,N
              DEPL(I)=DEPL(I)+R*SBAR(I,JP)
120       CONTINUE
200   CONTINUE
      IF (LHOOK) CALL DR_HOOK('DDR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE DDR
!
!-----------------------------------------------------------------------
!
