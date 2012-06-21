MODULE PRLE1AD_MOD
CONTAINS
SUBROUTINE PRLE1AD(KM,PLEPO)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM
USE TPM_DISTR
USE TPM_FIELDS
USE TPM_GEOMETRY

#ifdef DOC

!**** *PRLE1AD* - Prepare Legendre polonomials for inverse tranform.

!     Purpose.
!     --------
!        Copy the the Legendre polynomials for
!     specific zonal wavenumber KM to work arrays.

!**   Interface.
!     ----------
!        CALL PRLE1AD(...)

!        Explicit arguments :  KM - zonal wavenumber
!        -------------------   PLEPO - Legendre polonomial for zonal
!                                      wavenumber KM

!        Implicit arguments : 
!        --------------------

!     Method.
!     -------

!     Externals.   None
!     ----------   

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From PRLE1AD in IFS CY22R1
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN)  :: KM
REAL(KIND=JPRB),    INTENT(OUT) :: PLEPO(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) ::  ISL, JGL, JN


!     ------------------------------------------------------------------

!*       1.       COPY LEGENDRE POLONOMIALS.
!                 --------------------------

ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
DO JN=1,R%NSMAX-KM+2
  DO JGL=ISL,R%NDGNH
    PLEPO(JGL,JN) = F%RPNM(JGL,D%NPMS(KM)+JN)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE PRLE1AD
END MODULE PRLE1AD_MOD
