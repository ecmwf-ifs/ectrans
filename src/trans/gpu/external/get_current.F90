SUBROUTINE GET_CURRENT(KRESOL,LDLAM)

!**** *GET_CURRENT* - Extract current information from the transform package

!     Purpose.
!     --------
!     Interface routine for extracting current information from the T.P.

!**   Interface.
!     ----------
!     CALL GET_CURRENT(...)

!     Explicit arguments : (all optional)
!     -------------------- 
!     KRESOL   - Current resolution
!     LDLAM    -  .T. if the corresponding resolution is LAM, .F. if it is global

!     Method.
!     -------

!     Externals.  None
!     ----------  

!     Author.
!     -------
!        Ryad El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original : 24-Aug-2012

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM

!ifndef INTERFACE

USE TPM_GEN
USE TPM_GEOMETRY

!endif INTERFACE

IMPLICIT NONE

INTEGER(KIND=JPIM)  ,OPTIONAL,INTENT(OUT)  :: KRESOL
LOGICAL             ,OPTIONAL,INTENT(OUT)  :: LDLAM

!ifndef INTERFACE

! Get current resolution
IF (PRESENT(KRESOL)) KRESOL= NCUR_RESOL
IF (PRESENT(LDLAM))  LDLAM = G%LAM 


!endif INTERFACE

END SUBROUTINE GET_CURRENT
