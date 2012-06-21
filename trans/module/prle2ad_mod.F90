MODULE PRLE2AD_MOD
CONTAINS
SUBROUTINE PRLE2AD(KM,PLEPO)

!**** *PRLE2AD* - Prepare Legendre polonomials for direct transform.

!     Purpose.
!     --------
!        Copy the the Legendre polynomials for
!     specific zonal wavenumber M to work arrays. Multiply by
!     Gaussian weights for use in direct transforms.

!**   Interface.
!     ----------
!        CALL PRLE2AD(KM,PLEPO)

!        Explicit arguments :  KM - zonal wavenumber
!        -------------------   PLEPO - Legendre polonomial for zonal
!                                      wavenumber KM

!        Implicit arguments :  
!        --------------------

!     Method.
!     -------

!     Externals.   NONE
!     ----------  

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-27
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - NTMAX,NPMT instead of NSMAX,NPM
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_DISTR
USE TPM_FIELDS
USE TPM_GEOMETRY

IMPLICIT NONE

!     DUMMY ARGUMENTS
INTEGER_M, INTENT(IN)  :: KM
REAL_B,    INTENT(OUT) :: PLEPO(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: ISL, JGL, JN ,ISMAX, IDGNH, IPMT

!     ------------------------------------------------------------------

!*       1.    COPY LEGENDRE POLONOMIALS, MULTIPLY BY GAUSSIAN WEIGHTS.
!              --------------------------------------------------------

ISMAX = R%NSMAX
IDGNH = R%NDGNH
ISL   = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IPMT  = D%NPMT(KM)

DO JN=1,ISMAX-KM+2
  DO JGL=ISL,IDGNH
    PLEPO(JN,JGL)= F%RPNM(JGL,IPMT+JN)*F%RW(JGL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE PRLE2AD
END MODULE PRLE2AD_MOD
