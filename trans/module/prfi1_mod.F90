MODULE PRFI1_MOD
CONTAINS
SUBROUTINE PRFI1(KM,KF_UV,KF_SCALARS,PIA,PSPVOR,PSPDIV,PSPSCALAR,&
 & KFLDPTRUV,KFLDPTRSC)

#include "tsmbkind.h"

USE TPM_DISTR
USE TPM_TRANS

USE PRFI1B_MOD

#ifdef DOC

!**** *PRFI1* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1(KM,PIA,PSPVOR,PSPDIV,PSPSCALAR

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
!        Original : 00-02-01 From PRFI1 in IFS CY22R1

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KM,KF_UV,KF_SCALARS
REAL_B ,OPTIONAL, INTENT(IN)  :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(IN)  :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(IN)  :: PSPSCALAR(:,:)
REAL_B ,          INTENT(OUT) :: PIA(:,:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IDIV, IFIRST, ILAST, IVOR, JFLD


!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              ------------------------------------

IFIRST = 1
ILAST  = 4*KF_UV

!*       1.1      VORTICITY AND DIVERGENCE.

IF(KF_UV > 0)THEN
  IVOR = 1
  IDIV = 2*KF_UV+1
  CALL PRFI1B(KM,PIA(:,IVOR:IDIV-1),PSPVOR,KF_UV,KFLDPTRUV)
  CALL PRFI1B(KM,PIA(:,IDIV:ILAST) ,PSPDIV,KF_UV,KFLDPTRUV)
  ILAST = ILAST+4*KF_UV
ENDIF

!*       1.2    SCALAR VARIABLES.

IF(KF_SCALARS > 0)THEN
  IFIRST = ILAST+1
  ILAST  = IFIRST - 1 + 2*KF_SCALARS 
  CALL PRFI1B(KM,PIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE PRFI1
END MODULE PRFI1_MOD


