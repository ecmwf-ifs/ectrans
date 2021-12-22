SUBROUTINE DBFGSL (K_N,YD_D,K_M,K_NYS,K_JMIN,K_JMAX,YD_YBAR,YD_SBAR,P_RHO,P_SIZE)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND2  ,ONLY : JPRH

USE CONTROL_VECTORS_MOD

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: K_M 
INTEGER(KIND=JPIM)               :: K_N ! Argument NOT used
TYPE(CONTROL_VECTOR),INTENT(INOUT) :: YD_D 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_NYS 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JMIN 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JMAX 
TYPE(CONTROL_VECTOR),INTENT(IN)    :: YD_YBAR(K_M) 
TYPE(CONTROL_VECTOR),INTENT(IN)    :: YD_SBAR(K_M) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_RHO(K_M) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SIZE 
!----

!     compute the product H u, where
!     . H is the matrix that could be obtained from the m pairs
!       (ybar,sbar) using the inverse BFGS formula from the diagonal
!       matrix (size) * (Identity matrix),
!     . u is a vector of dimension n.

!     The vector d contains

!       u (on entry), H u (on output).

!     rho(m) is a working zone.

!     Input:

!----

! --- local variables

INTEGER(KIND=JPIM) :: J,JFIN,JP
REAL(KIND=JPRB) :: Z_R
REAL(KIND=JPRH) :: Z_S
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "dpseuclid.h"

! --- return if there is no pair (y,s) in ybar and sbar

IF (LHOOK) CALL DR_HOOK('DBFGSL',0,ZHOOK_HANDLE)
IF (K_NYS == 0 .AND. LHOOK) CALL DR_HOOK('DBFGSL',1,ZHOOK_HANDLE)
IF (K_NYS == 0) RETURN

! --- set jfin

JFIN=K_JMAX
IF (JFIN < K_JMIN) JFIN=K_JMAX+K_M

! --- backward sweep

DO J=JFIN,K_JMIN,-1
  JP=J
  IF (JP > K_M) JP=JP-K_M
  CALL DPSEUCLID (K_N,YD_D,YD_SBAR(JP),Z_S)
  P_RHO(JP)=Z_S
  YD_D%DATA = YD_D%DATA - Z_S*YD_YBAR(JP)%DATA
ENDDO

! --  preconditioning

Z_R=P_SIZE
YD_D%DATA = YD_D%DATA*Z_R

! --- forward sweep

DO J=K_JMIN,JFIN
  JP=J
  IF (JP > K_M) JP=JP-K_M
  CALL DPSEUCLID (K_N,YD_D,YD_YBAR(JP),Z_S)
  Z_R=P_RHO(JP)-Z_S
  YD_D%DATA = YD_D%DATA + Z_R*YD_SBAR(JP)%DATA
ENDDO
IF (LHOOK) CALL DR_HOOK('DBFGSL',1,ZHOOK_HANDLE)
END SUBROUTINE DBFGSL
