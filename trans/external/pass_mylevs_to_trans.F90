SUBROUTINE PASS_MYLEVS_TO_TRANS(KMYLEVS)

!**** Pass MYLEVS array to transform package, necessary for call to FSPGLH
!     KMYLEVS         - list of levels handled by this processor/PE, necessary for call to FSPGLH

!     -------
!     Original version:    O. Marsden *ECMWF*  Aug 2016
!     -------

USE PARKIND1  ,ONLY : JPIM,JPRB 

!ifndef INTERFACE

USE TPM_DISTR       ,ONLY : D
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!endif INTERFACE

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM) , INTENT(IN) :: KMYLEVS(:)

!ifndef INTERFACE

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "user_clock.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PASS_MYLEVS_TO_TRANS',0,ZHOOK_HANDLE)


  IF (.NOT. ALLOCATED(D%NMYLEVS)) THEN
    ! Allocate and fill MYLEVS array for D distr_type element
    ALLOCATE(D%NMYLEVS(SIZE(KMYLEVS)))
    D%NMYLEVS = KMYLEVS
  ELSE
    CALL ABORT_TRANS('PASS_MYLEVS_TO_TRANS: Attempting to allocate D%NMYLEVS more than once - regretfully aborting')
  END IF

IF (LHOOK) CALL DR_HOOK('PASS_MYLEVS_TO_TRANS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE PASS_MYLEVS_TO_TRANS


