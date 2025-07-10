SUBROUTINE ETRANS_RELEASE(KRESOL)

!**** *ETRANS_RELEASE* - release a spectral resolution

!     Purpose.
!     --------
!      Release all arrays related to a given resolution tag

!**   Interface.
!     ----------
!     CALL ETRANS_RELEASE

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

USE EDEALLOC_RESOL_MOD   ,ONLY : EDEALLOC_RESOL
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KRESOL

!endif INTERFACE

!     ------------------------------------------------------------------

CALL EDEALLOC_RESOL(KRESOL)

!     ------------------------------------------------------------------

END SUBROUTINE ETRANS_RELEASE
