SUBROUTINE TRANS_RELEASE(KRESOL)

!**** *TRANS_RELEASE* - release a spectral resolution

!     Purpose.
!     --------
!      Release all arrays related to a given resolution tag

!**   Interface.
!     ----------
!     CALL TRANS_RELEASE

!     Explicit arguments : KRESOL : resolution tag
!     --------------------

!     Method.
!     -------

!     Externals.  None
!     ----------

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 09-Jul-2013

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM

!ifndef INTERFACE

USE DEALLOC_RESOL_MOD   ,ONLY : DEALLOC_RESOL
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KRESOL

!endif INTERFACE

!     ------------------------------------------------------------------

CALL DEALLOC_RESOL(KRESOL)

!     ------------------------------------------------------------------

END SUBROUTINE TRANS_RELEASE
