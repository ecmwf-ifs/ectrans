INTERFACE
SUBROUTINE TRANS_END(CDMODE)

!**** *TRANS_END* - Terminate transform package 

!     Purpose.
!     --------
!     Terminate transform package. Release all allocated arrays.

!**   Interface.
!     ----------
!     CALL TRANS_END

!     Explicit arguments : None
!     -------------------- 

!     Method.
!     -------

!     Externals.  None
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03


!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
CHARACTER*5, OPTIONAL, INTENT(IN) :: CDMODE


END SUBROUTINE TRANS_END
END INTERFACE
