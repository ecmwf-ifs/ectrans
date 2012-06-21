MODULE LEDIR_MOD
CONTAINS
SUBROUTINE LEDIR(KM,KFC,KLED2,PAIA,PSIA,POA1,PLEPO)

!**** *LEDIR* - Direct Legendre transform.

!     Purpose.
!     --------
!        Direct Legendre tranform of state variables.

!**   Interface.
!     ----------
!        CALL LEDIR(...)

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

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_TRANS

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M, INTENT(IN)  :: KM,KFC,KLED2


REAL_B ,   INTENT(IN)  :: PSIA(:,:),   PAIA(:,:)
REAL_B,    INTENT(IN)  :: PLEPO(:,:)
REAL_B :: POA1(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IA, IDGLU, IFC, ILA, ILS, IS, ISKIP, ISL


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
ELSE
  ISKIP = 1
  IFC   = KFC
ENDIF


IF (IFC > 0) THEN
!*       1.2      ANTISYMMETRIC PART.

  IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
  CALL MXMAOP(PLEPO(IA,ISL),2,R%NLED3,PAIA(1,ISL),KLED2,ISKIP &
      &,POA1(IA,1),2,R%NLED4*ISKIP,ILA,IDGLU,IFC)

!*       1.3      SYMMETRIC PART.

  CALL MXMAOP(PLEPO(IS,ISL),2,R%NLED3,PSIA(1,ISL),KLED2,ISKIP &
      &,POA1(IS,1),2,R%NLED4*ISKIP,ILS,IDGLU,IFC)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
