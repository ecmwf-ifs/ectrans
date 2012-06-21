MODULE LEDIRAD_MOD
CONTAINS
SUBROUTINE LEDIRAD(KM,KFC,KLED2,PAIA,PSIA,POA1,PLEPO)

!**** *LEDIRAD* - Direct Legendre transform.

!     Purpose.
!     --------
!        Direct Legendre tranform of state variables.

!**   Interface.
!     ----------
!        CALL LEDIRAD(...)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFC - number of field to transform
!                              PAIA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSIA - symmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              POA1 -  spectral
!                              fields for zonal wavenumber KM
!                              PLEPO - Legendre polonomials

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - matrix multiply
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - NTMAX instead of NSMAX
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_TRANS

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KLED2

REAL(KIND=JPRB),    INTENT(OUT)  :: PSIA(:,:),   PAIA(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PLEPO(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: POA1(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IA, IDGLU, IFC, ILA, ILS, IS, ISKIP, ISL, JGL, J1


!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

IA  = 1+MOD(R%NTMAX-KM+2,2)
IS  = 1+MOD(R%NTMAX-KM+1,2)
ILA = (R%NTMAX-KM+2)/2
ILS = (R%NTMAX-KM+3)/2
ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

IF(KM == 0)THEN
  ISKIP = 2
  IFC   = KFC/2
  DO JGL=ISL,R%NDGNH
    DO J1=2,KFC,2
      PSIA(J1,JGL)=0.0_JPRB
      PAIA(J1,JGL)=0.0_JPRB
    ENDDO
  ENDDO
ELSE
  ISKIP = 1
  IFC   = KFC
ENDIF


IF (IFC > 0) THEN
!*       1.2      ANTISYMMETRIC PART.

  IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
  CALL MXMAOP(PLEPO(IA,ISL),R%NLED3,2,POA1(IA,1),2,R%NLED4*ISKIP, &
      & PAIA(1,ISL),KLED2,ISKIP,IDGLU,ILA,IFC)

!*       1.3      SYMMETRIC PART.

  CALL MXMAOP(PLEPO(IS,ISL),R%NLED3,2,POA1(IS,1),2,R%NLED4*ISKIP, &
      & PSIA(1,ISL),KLED2,ISKIP,IDGLU,ILS,IFC)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LEDIRAD
END MODULE LEDIRAD_MOD
