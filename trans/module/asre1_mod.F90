MODULE ASRE1_MOD
CONTAINS
SUBROUTINE ASRE1(KM,KMLOC,PAOA1,PSOA1)

#include "tsmbkind.h"

USE TPM_TRANS
 
USE ASRE1B_MOD


!**** *ASRE1* - Recombine antisymmetric and symmetric parts

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!       *CALL* *ASRE1(...)

!        Explicit arguments : 
!        --------------------  
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (basic
!                              variables and N-S derivatives)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (basic
!                              variables and N-S derivatives)

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.   ASRE1B - basic recombination routine
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From ASRE1 in IFS CY22R1

!     ------------------------------------------------------------------


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M , INTENT(IN) :: KM
INTEGER_M , INTENT(IN) :: KMLOC

REAL_B    , INTENT(IN) :: PSOA1(:,:),       PAOA1(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IFLDS


!     ------------------------------------------------------------------

IFLDS = NF_OUT_LT

CALL ASRE1B(IFLDS,KM,KMLOC,PAOA1,PSOA1)

!     ------------------------------------------------------------------

END SUBROUTINE ASRE1
END MODULE ASRE1_MOD
