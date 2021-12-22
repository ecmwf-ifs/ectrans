      SUBROUTINE DD (CTONB,CTCAB,N,SSCALE,NM,DEPL,AUX,JMIN,JMAX,        &
     &               PRECOS,DIAG,ALPHA,YBAR,SBAR)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!----
!
!     calcule le produit H.g ou
!         . H est une matrice construite par la formule de bfgs inverse
!           a nm memoires a partir de la matrice diagonale diag
!           dans un espace hilbertien dont le produit scalaire
!           est donne par prosca
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
!     izs(1),rzs(1),dzs(1) sont des zones de travail pour prosca
!
!     Modfications:
!     LF. Meunier       05-03-2015: Add some OpenMP
!----
      USE CONTROL_VECTORS_MOD
      IMPLICIT NONE
!
!         arguments
!
      LOGICAL :: SSCALE
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: NM
      INTEGER(KIND=JPIM) :: JMIN
      INTEGER(KIND=JPIM) :: JMAX
      REAL(KIND=JPRB) :: PRECOS
      REAL(KIND=JPRB) :: ALPHA(NM)
      TYPE (CONTROL_VECTOR) :: DEPL
      TYPE (CONTROL_VECTOR) :: DIAG
      TYPE (CONTROL_VECTOR) :: AUX
      TYPE (CONTROL_VECTOR),DIMENSION(NM) :: YBAR
      TYPE (CONTROL_VECTOR),DIMENSION(NM) :: SBAR
      EXTERNAL :: CTONB,CTCAB
!
!         variables locales
!
      INTEGER(KIND=JPIM) :: JFIN,I,J,JP
      REAL(KIND=JPRB) :: R,PS

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('DD',0,ZHOOK_HANDLE)
      CALL CTONB (N,DEPL,AUX)

      JFIN=JMAX
      IF (JFIN.LT.JMIN) JFIN=JMAX+NM
!
!         phase de descente
!
      DO J=JFIN,JMIN,-1
          JP=J
          IF (JP.GT.NM) JP=JP-NM
          PS = DOT_PRODUCT (AUX,SBAR(JP))
          R=PS
          ALPHA(JP)=R
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(I)
          DO I=1, AUX%NSIZEL
            AUX%DATA(I) = AUX%DATA(I) -R * YBAR(JP)%DATA(I)
          ENDDO
!$OMP END PARALLEL DO
      ENDDO
!
!         preconditionnement
!


      IF (SSCALE) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(I)
        DO I=1, AUX%NSIZEl
          AUX%DATA(I) = AUX%DATA(I) * PRECOS
        ENDDO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(I)
        DO I=1, AUX%NSIZEL
          AUX%DATA(I) = AUX%DATA(I) * DIAG%DATA(I)
        ENDDO
!$OMP END PARALLEL DO
      ENDIF



!
!         remontee
!
      DO J=JMIN,JFIN
          JP=J
          IF (JP.GT.NM) JP=JP-NM
          PS = DOT_PRODUCT (AUX,YBAR(JP))
          R=ALPHA(JP)-PS
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(I)
          DO I=1, AUX%NSIZEL 
            AUX%DATA(I) = AUX%DATA(I) + R * SBAR(JP)%DATA(I)
          ENDDO
!$OMP END PARALLEL DO
      ENDDO

      CALL CTCAB (N,AUX,DEPL)
      IF (LHOOK) CALL DR_HOOK('DD',1,ZHOOK_HANDLE)
      END SUBROUTINE DD
!
!-----------------------------------------------------------------------
!
