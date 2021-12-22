SUBROUTINE MULTVDV(PVEC,PDIA,PROD)

!     Purpose.
!     --------
!       Compute A=VDV^t

!        Explicit arguments :
!        --------------------
!         PVEC : matrix V           (in)
!         PDIA : diagonal D         (in)
!         PROD : matrix product A   (out)

!     Author.
!     -------
!       Yannick Tremolet

!     Modifications.
!     --------------
!       Original  : 12-May-2004
! ------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK , ONLY: LHOOK, DR_HOOK

! ------------------------------------------------------------------

IMPLICIT NONE
REAL(KIND=JPRB), INTENT(IN)  :: PVEC(:,:)
REAL(KIND=JPRB), INTENT(IN)  :: PDIA(:)
REAL(KIND=JPRB), INTENT(OUT) :: PROD(:,:)

! ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILEN,JI,JJ,JK
REAL(KIND=JPRB), ALLOCATABLE :: ZZ(:,:)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! ------------------------------------------------------------------

#include "abor1.intfb.h"

! ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MULTVDV',0,ZHOOK_HANDLE)
! ------------------------------------------------------------------

ILEN=SIZE(PDIA)
IF (SIZE(PVEC,1)/=ILEN) CALL ABOR1('MULVDV: Error size PVEC 1')
IF (SIZE(PVEC,2)/=ILEN) CALL ABOR1('MULVDV: Error size PVEC 2')
IF (SIZE(PROD,1)/=ILEN) CALL ABOR1('MULVDV: Error size PROD 1')
IF (SIZE(PROD,2)/=ILEN) CALL ABOR1('MULVDV: Error size PROD 2')

ALLOCATE(ZZ(ILEN,ILEN))
ZZ(:,:)=0.0_JPRB

DO JJ=1,ILEN
  DO JI=1,ILEN
    ZZ(JI,JJ)=ZZ(JI,JJ)+PVEC(JI,JJ)*PDIA(JJ)
  ENDDO
ENDDO

PROD(:,:)=0.0_JPRB

! Diagonal and above
DO JK=1,ILEN
  DO JJ=1,ILEN
    DO JI=1,JJ
      PROD(JI,JJ)=PROD(JI,JJ)+ZZ(JI,JK)*PVEC(JJ,JK)
    ENDDO
  ENDDO
ENDDO

! Below diagonal using symmetry
DO JJ=1,ILEN-1
  DO JI=JJ+1,ILEN
    PROD(JI,JJ)=PROD(JJ,JI)
  ENDDO
ENDDO

DEALLOCATE(ZZ)

! ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MULTVDV',1,ZHOOK_HANDLE)
END SUBROUTINE MULTVDV
