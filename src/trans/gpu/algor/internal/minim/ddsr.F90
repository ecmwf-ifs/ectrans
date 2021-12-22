      SUBROUTINE DDSR (YDDIM,PROSCAR,CTONBR,CTCABR,N,SSCALE,NM,DEPL,AUX,&
     &                 JMIN,JMAX,PRECOS,DIAG,ALPHA,YBAR,SBAR)
!AUTOPROMOTE

      USE YOMDIM  , ONLY : TDIM
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!----
!
!     This subroutine has the same role as ddr (computation of the
!     product H.g). It supposes however that the (y,s) pairs are not
!     stored in core memory, but on a devise chosen by the user.
!     The access to this devise is performed via the subroutine ystbl.
!
!----
!
!         arguments
!
      IMPLICIT NONE
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
      REAL(KIND=JPRB) :: YBAR(N)
      REAL(KIND=JPRB) :: SBAR(N)
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
      IF (LHOOK) CALL DR_HOOK('DDSR',0,ZHOOK_HANDLE)
      JFIN=JMAX
      IF (JFIN.LT.JMIN) JFIN=JMAX+NM
!
!         phase de descente
!
      DO 100 J=JFIN,JMIN,-1
          JP=J
          IF (JP.GT.NM) JP=JP-NM
          CALL YSTBLR (.FALSE.,YBAR,SBAR,N,JP)
          CALL PROSCAR (YDDIM,N,DEPL,SBAR,PS)
          R=PS
          ALPHA(JP)=R
          DO 20 I=1,N
              DEPL(I)=DEPL(I)-R*YBAR(I)
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
          CALL YSTBLR (.FALSE.,YBAR,SBAR,N,JP)
          CALL PROSCAR (YDDIM,N,DEPL,YBAR(1),PS)
          R=ALPHA(JP)-PS
          DO 120 I=1,N
              DEPL(I)=DEPL(I)+R*SBAR(I)
120       CONTINUE
200   CONTINUE
      IF (LHOOK) CALL DR_HOOK('DDSR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE DDSR
!
!-----------------------------------------------------------------------
!
