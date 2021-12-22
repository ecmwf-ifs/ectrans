      SUBROUTINE DDS (CTONB,CTCAB,N,SSCALE,NM,DEPL,AUX,JMIN,JMAX,       &
     &                PRECOS,DIAG,ALPHA,YBAR,SBAR)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!----
!
!     This subroutine has the same role as dd (computation of the
!     product H.g). It supposes however that the (y,s) pairs are not
!     stored in core memory, but on a devise chosen by the user.
!     The access to this devise is performed via the subroutine ystbl.
!
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
      TYPE (CONTROL_VECTOR) :: YBAR
      TYPE (CONTROL_VECTOR) :: SBAR
      EXTERNAL :: CTONB,CTCAB
!
!         variables locales
!
      INTEGER(KIND=JPIM) :: JFIN,I,J,JP
      REAL(KIND=JPRB) :: R,PS

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('DDS',0,ZHOOK_HANDLE)
      CALL CTONB (N,DEPL,AUX)

      JFIN=JMAX
      IF (JFIN.LT.JMIN) JFIN=JMAX+NM
!
!         phase de descente
!
      DO J=JFIN,JMIN,-1
          JP=J
          IF (JP.GT.NM) JP=JP-NM
          STOP 'ystbl:dead code'
!         call ystbl (.false.,ybar,sbar,n,jp)
          PS = DOT_PRODUCT (AUX,SBAR)
          R=PS
          ALPHA(JP)=R
          AUX%DATA = AUX%DATA - R * YBAR%DATA
      ENDDO
!
!         preconditionnement
!
      IF (SSCALE) THEN
          AUX%DATA = AUX%DATA * PRECOS
      ELSE
          AUX%DATA = AUX%DATA * DIAG%DATA
      ENDIF
!
!         remontee
!
      DO J=JMIN,JFIN
          JP=J
          IF (JP.GT.NM) JP=JP-NM
          STOP 'ystbl:dead code'
!         call ystbl (.false.,ybar,sbar,n,jp)
          PS = DOT_PRODUCT (AUX,YBAR)
          R=ALPHA(JP)-PS
          AUX%DATA = AUX%DATA + R * SBAR%DATA
      ENDDO

      CALL CTCAB (N,AUX,DEPL)
      IF (LHOOK) CALL DR_HOOK('DDS',1,ZHOOK_HANDLE)
      ENDSUBROUTINE DDS
!
!-----------------------------------------------------------------------
!
