MODULE LDFOU2_MOD
CONTAINS
SUBROUTINE LDFOU2(KM,KF_UV,PAIA,PSIA)

!**** *LDFOU2* - Division by a*cos(theta) of u and v

!     Purpose.
!     --------
!        In Fourier space divide u and v by  a*cos(theta).

!**   Interface.
!     ----------
!        CALL LDFOU2(KM,PAIA,PSIA)

!        Explicit arguments : 
!        --------------------  KM - zonal wavenumber
!                              PAIA - antisymmetric fourier fields
!                              PSIA - symmetric fourierfields

!        Implicit arguments :  RACTHE - 1./(a*cos(theta))
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 91-07-01
!        Modified : 94-04-06 R. El Khatib - Full-POS configuration 'P'
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                    instead of u,v->vor,div
!        MPP Group: 95-10-01 Message Passing option added
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_FIELDS

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M, INTENT(IN) :: KM,KF_UV

REAL_B ,INTENT(INOUT) :: PSIA(:,:),   PAIA(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: J, JGL ,IFLD ,ISL


!     ------------------------------------------------------------------

!*       1.    DIVIDE U V BY A*COS(THETA)
!              --------------------------

ISL  = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IFLD = 4*KF_UV

!*       1.1      U AND V 

DO JGL=ISL,R%NDGNH
  DO J=1,IFLD
    PAIA(J,JGL) = PAIA(J,JGL)*F%RACTHE(JGL)
    PSIA(J,JGL) = PSIA(J,JGL)*F%RACTHE(JGL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE LDFOU2
END MODULE LDFOU2_MOD
