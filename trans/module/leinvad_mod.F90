MODULE LEINVAD_MOD
CONTAINS
SUBROUTINE LEINVAD(KM,KFC,PIA,PAOA1,PSOA1,PLEPO)

!**** *LEINVAD* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL LEINVAD(...)

!        Explicit arguments :  KM - zonal wavenumber (input-c)
!        --------------------  KFC - number of fields to tranform (input-c)
!                              PIA - spectral fields
!                              for zonal wavenumber KM (input)
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PLEPO - Legendre polonomials for zonal
!                              wavenumber KM (input-c)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - calls SGEMVX (matrix multiply)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LEINVAD in IFS CY22R1

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_TRANS

IMPLICIT NONE

INTEGER_M, INTENT(IN)    :: KM,KFC
REAL_B,    INTENT(OUT)   :: PIA(:,:)
REAL_B,    INTENT(IN)    :: PLEPO(:,:)
REAL_B,    INTENT(INOUT) :: PSOA1(:,:)
REAL_B,    INTENT(INOUT) :: PAOA1(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IA, IDGLU, IFC, ILA, ILS, IS, ISKIP, ISL, J1, JGL,IOAD1


!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

IA  = 1+MOD(R%NSMAX-KM+2,2)
IS  = 1+MOD(R%NSMAX-KM+1,2)
ILA = (R%NSMAX-KM+2)/2
ILS = (R%NSMAX-KM+3)/2
ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IOAD1 = 2*NF_OUT_LT

IF(KM == 0)THEN
  ISKIP = 2
  IFC = KFC/2
ELSE
  ISKIP = 1
  IFC = KFC
ENDIF

!*       1.2      ANTISYMMETRIC PART.

IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
CALL MXMAOP(PLEPO(ISL,IA),2*R%NLEI3,1,PAOA1(1,ISL),IOAD1,&
            &ISKIP,PIA(1+IA,1),2,R%NLEI1*ISKIP,&
            &ILA,IDGLU,IFC)

!*       1.3      SYMMETRIC PART.

CALL MXMAOP(PLEPO(ISL,IS),2*R%NLEI3,1,PSOA1(1,ISL),IOAD1,&
            &ISKIP,PIA(1+IS,1),2,R%NLEI1*ISKIP,&
            &ILS,IDGLU,IFC)

!     ------------------------------------------------------------------


END SUBROUTINE LEINVAD
END MODULE LEINVAD_MOD
