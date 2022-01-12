INTERFACE
SUBROUTINE TRANS_PNM(KRESOL,KM,PRPNM,LDTRANSPOSE,LDCHEAP)

!**** *TRANS_PNM* - Compute Legendre polynomials for a given wavenember

!     Purpose.
!     --------
!     Interface routine for computing Legendre polynomials for a given wavenember

!**   Interface.
!     ----------
!     CALL TRANS_PNM(...)

!     Explicit arguments : All arguments are optional.
!     --------------------
!     KRESOL   - resolution tag for which info is required ,default is the
!                first defined resulution (input)
!     KM       - wave number
!     PRPNM    - Legendre polynomials
!     LDTRANSPOSE - Legendre polynomials array is transposed
!     LDCHEAP   - cheapest but less accurate computation

!     Method.
!     -------

!     Externals.  SET_RESOL - set resolution
!     ----------

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 22-Jan-2016

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KM
REAL(KIND=JPRB)    ,INTENT(OUT) :: PRPNM(:,:)
LOGICAL, OPTIONAL, INTENT(IN) :: LDTRANSPOSE
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHEAP

END SUBROUTINE TRANS_PNM
END INTERFACE
