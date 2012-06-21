MODULE UPDSPAD_MOD
CONTAINS
SUBROUTINE UPDSPAD(KM,KF_UV,KF_SCALARS,POA1,POA2, &
 & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC)

!**** *UPDSPAD* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update the spectral arrays for a fixed zonal wave-number
!        from values in POA1 and POA2.

!**   Interface.
!     ----------
!        CALL UPDSPAD(...)

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

!     Externals.  UPDSPADB - basic transfer routine
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
!        Modified : 94-08-02 R. El Khatib - interface to UPDSPADB
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR

USE UPDSPBAD_MOD

IMPLICIT NONE


!     DUMMY INTEGER SCALARS

INTEGER_M, INTENT(IN)  :: KM,KF_UV,KF_SCALARS

REAL_B , INTENT(OUT)  :: POA1(:,:)
REAL_B , INTENT(OUT)  :: POA2(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
INTEGER_M,OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IVORS, IVORE, IDIVS, IDIVE, IST ,IEND, JN, ISE


!     ------------------------------------------------------------------

!*       1.    UPDATE FIELDS
!              -------------


!*       1.1      VORTICITY AND DIVERGENCE.

IF (KF_UV > 0) THEN
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV
  IF (KM == 0) THEN
    DO JN=0,R%NSMAX
      ISE = 1+JN*2+1
      PSPDIV(:,ISE) = _ZERO_
      PSPVOR(:,ISE) = _ZERO_
    ENDDO
  ENDIF
  CALL UPDSPBAD(KM,KF_UV,POA2(:,IVORS:IVORE),PSPVOR,KFLDPTRUV)
  CALL UPDSPBAD(KM,KF_UV,POA2(:,IDIVS:IDIVE),PSPDIV,KFLDPTRUV)
ENDIF

!*       1.2   SCALARS

IF (KF_SCALARS > 0) THEN
  IST = 1
  IEND = 2*KF_SCALARS
  IF(KF_UV > 0) THEN
    IST = IST+4*KF_UV
    IEND = IEND+4*KF_UV
  ENDIF
  IF (KM == 0) THEN
    DO JN=0,R%NSMAX
      ISE = 1+JN*2+1
      PSPSCALAR(:,ISE) = _ZERO_
    ENDDO
  ENDIF
  CALL UPDSPBAD(KM,KF_SCALARS,POA1(:,IST:IEND),PSPSCALAR,KFLDPTRSC)

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UPDSPAD
END MODULE UPDSPAD_MOD
