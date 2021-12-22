      SUBROUTINE YSTBL_1DV (STORE,YBAR,SBAR,N,J)
!----
!
!     This subroutine should store (if store = .true.) or restore
!     (if store = .false.) a pair (ybar,sbar) at or from position
!     j in memory. Be sure to have 1 <= j <= m, where m in the number
!     of updates specified by subroutine mupdts.
!
!     The subroutine is used only when the (y,s) pairs are not
!     stored in core memory in the arrays ybar(.,.) and sbar(.,.).
!     In this case, the subroutine has to be written by the user.
!
!----
!
!         arguments
!
      USE PARKIND1, ONLY : JPIM, JPRB
!
      IMPLICIT NONE
!
      LOGICAL :: STORE
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: J
      REAL(KIND=JPRB) :: YBAR(N)
      REAL(KIND=JPRB) :: SBAR(N)
!
!RJ: does not feel right
      CALL ABOR1("where I in ystblr_1dv? zRJ") !RJ-debug
      RETURN
      ENDSUBROUTINE YSTBL_1DV
