MODULE UPDSP_MOD
CONTAINS
SUBROUTINE UPDSP(KM,POA1,POA2,PSPVOR,PSPDIV,PSPSCALAR)

!**** *UPDSP* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update the spectral arrays for a fixed zonal wave-number
!        from values in POA1 and POA2.

!**   Interface.
!     ----------
!        CALL UPDSP(...)

!        Explicit arguments : 
!        -------------------- 
!        KM - zonal wave-number
!        POA1 - spectral fields for zonal wavenumber KM (basic var.)
!        POA2 - spectral fields for zonal wavenumber KM (vor. div.)
!        PSPVOR - spectral vorticity
!        PSPDIV - spectral divergence
!        PSPSCALAR - spectral scalar variables

!        Implicit arguments :  
!        --------------------

!     Method.
!     -------

!     Externals.  UPDSPB - basic transfer routine
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-02-02
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 94-08-02 R. El Khatib - interface to UPDSPB
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_TRANS
USE TPM_DISTR

USE UPDSPB_MOD

IMPLICIT NONE


!     DUMMY INTEGER SCALARS

INTEGER_M, INTENT(IN)  :: KM

REAL_B , INTENT(IN)  :: POA1(:,:)
REAL_B , INTENT(IN)  :: POA2(:,:)
REAL_B  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IVORS, IVORE, IDIVS, IDIVE, IST ,IEND


!     ------------------------------------------------------------------

!*       1.    UPDATE FIELDS
!              -------------


!*       1.1      VORTICITY AND DIVERGENCE.

IF (NF_UV > 0) THEN
  IVORS = 1
  IVORE = 2*NF_UV
  IDIVS = 2*NF_UV+1
  IDIVE = 4*NF_UV
  CALL UPDSPB(KM,NF_UV,POA2(:,IVORS:IVORE),PSPVOR)
  CALL UPDSPB(KM,NF_UV,POA2(:,IDIVS:IDIVE),PSPDIV)
ENDIF

!*       1.2   SCALARS

IF (NF_SCALARS > 0) THEN
  IST = 1
  IEND = 2*NF_SCALARS
  IF(NF_UV > 0) THEN
    IST = IST+4*NF_UV
    IEND = IEND+4*NF_UV
  ENDIF
  CALL UPDSPB(KM,NF_SCALARS,POA1(:,IST:IEND),PSPSCALAR)

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UPDSP
END MODULE UPDSP_MOD
