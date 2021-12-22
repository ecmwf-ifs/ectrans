      SUBROUTINE MUPDTSR (SSCALE,INMEMO,N,M,NRZ)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
!         arguments
!
      IMPLICIT NONE
!
      INTEGER(KIND=JPIM) :: NRZ
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: M
      LOGICAL :: SSCALE
      LOGICAL :: INMEMO
!----
!
!     This routine has to return:
!       m:      the numer of updates to form the approximate Hessien H,
!       inmemo = .true.  if the vectors y and s used to form H are
!                        stored in core memory,
!                .false. otherwise,
!     When inmemo=.false., the routine `ystbl', which stores and
!     restores (y,s) pairs, has to be rewritten.
!
!----
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('MUPDTSR',0,ZHOOK_HANDLE)
      M=(NRZ-4*N)/(2*N+1)
      INMEMO=.TRUE.
      IF (LHOOK) CALL DR_HOOK('MUPDTSR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE MUPDTSR
