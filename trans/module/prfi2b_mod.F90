MODULE PRFI2B_MOD
CONTAINS
SUBROUTINE PRFI2B(KFIELD,KM,KMLOC,PAIA,PSIA)

!**** *PRFI2B* - Prepare input work arrays for direct transform

!     Purpose.
!     --------
!        To extract the Fourier fields for a specific zonal wavenumber
!        and put them in an order suitable for the direct Legendre
!        tranforms, i.e. split into symmetric and anti-symmetric part.

!**   Interface.
!     ----------
!     *CALL* *PRFI2B(..)

!        Explicit arguments : 
!        -------------------   KFIELD - number of fields
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM

!        Implicit arguments :  FOUBUF in TPM_TRANS
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 90-07-01
!        MPP Group: 95-10-01 Support for Distributed Memory version
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM
USE TPM_TRANS
USE TPM_GEOMETRY
USE TPM_DISTR

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD,KM,KMLOC
REAL(KIND=JPRB)  , INTENT(OUT) :: PSIA(:,:),   PAIA(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IGLS,  ISL, ISTAN, ISTAS, JF, JGL


!     ------------------------------------------------------------------

!*       1.    EXTRACT SYM./ANTISYM. FIELDS FROM FOURIER ARRAY.
!              ------------------------------------------------

ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

DO JGL=ISL,R%NDGNH
  IGLS = R%NDGL+1-JGL
  ISTAN = (D%NSTAGT1B(D%NPROCL(JGL) )+D%NPNTGTB1(KMLOC,JGL ))*2*KFIELD
  ISTAS = (D%NSTAGT1B(D%NPROCL(IGLS))+D%NPNTGTB1(KMLOC,IGLS))*2*KFIELD
!DIR$ IVDEP
!OCL    NOVREC
  DO JF=1,KFIELD*2
    PSIA(JF,JGL) = FOUBUF(ISTAN+JF)+FOUBUF(ISTAS+JF)
    PAIA(JF,JGL) = FOUBUF(ISTAN+JF)-FOUBUF(ISTAS+JF)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE PRFI2B
END MODULE PRFI2B_MOD
