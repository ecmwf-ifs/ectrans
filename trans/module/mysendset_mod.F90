MODULE MYSENDSET_MOD
CONTAINS
FUNCTION MYSENDSET(KSETS,KMYSET,KSET)

#ifdef DOC

!**** *MYSENDSET* RETURNS SET NUMBER TO SEND TO

!     Purpose.
!     --------
!       

!**   Interface.
!     ----------
!        ISENDSET = MYSENDSET(KSETS,KMYSET,KSET)

!        Explicit arguments :  
!        --------------------
!                  input:   KSETS

!        Implicit arguments :  NONE
!        --------------------
!     Method.
!     -------

!     Externals.
!     ----------
!         NONE

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-03

!     ------------------------------------------------------------------
#endif

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE ABORT_TRANS_MOD

IMPLICIT NONE
INTEGER(KIND=JPIM) :: MYSENDSET
INTEGER(KIND=JPIM),INTENT(IN)  :: KSETS,KMYSET,KSET


!     ------------------------------------------------------------------

!*       1.    Check input argument for validity 
!              ---------------------------------

IF(KSETS < 1 .OR. KMYSET > KSETS .OR. KSET > KSETS-1) THEN

  CALL ABORT_TRANS(' MYSENDSET: INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  MYSENDSET = MOD(KMYSET+KSET-1,KSETS)+1

ENDIF

END FUNCTION MYSENDSET
END MODULE MYSENDSET_MOD
