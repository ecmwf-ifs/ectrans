MODULE PRFI1AD_MOD
CONTAINS
SUBROUTINE PRFI1AD(KM,PIA,PSPVOR,PSPDIV,PSPSCALAR)

#include "tsmbkind.h"

USE TPM_DISTR
USE TPM_TRANS

USE PRFI1BAD_MOD

#ifdef DOC

!**** *PRFI1AD* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1AD(KM,PIA,PSPVOR,PSPDIV,PSPSCALAR

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              PSPVOR - vorticity
!                              PSPDIV - divergence
!                              PSPSCALAR - scalar variables

!        Implicit arguments :  None.
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
!        Original : 00-02-01 From PRFI1AD in IFS CY22R1

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KM
REAL_B ,OPTIONAL, INTENT(INOUT)  :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(INOUT)  :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(INOUT)  :: PSPSCALAR(:,:)
REAL_B ,          INTENT(IN)     :: PIA(:,:)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IDIV, IFIRST, ILAST, IVOR


!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              ------------------------------------

IFIRST = 1
ILAST  = 4*NF_UV

!*       1.1      VORTICITY AND DIVERGENCE.

IF(NF_UV > 0)THEN
  IVOR = 1
  IDIV = 2*NF_UV+1
  CALL PRFI1BAD(KM,PIA(:,IVOR:IDIV-1),PSPVOR,NF_UV)
  CALL PRFI1BAD(KM,PIA(:,IDIV:ILAST) ,PSPDIV,NF_UV)
  ILAST = ILAST+4*NF_UV
ENDIF

!*       1.2    SCALAR VARIABLES.

IF(NF_SCALARS > 0)THEN
  IFIRST = ILAST+1
  ILAST  = IFIRST - 1 + 2*NF_SCALARS 
  CALL PRFI1BAD(KM,PIA(:,IFIRST:ILAST),PSPSCALAR(:,:),NF_SCALARS)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE PRFI1AD
END MODULE PRFI1AD_MOD


