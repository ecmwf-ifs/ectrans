MODULE LTINV_MOD
CONTAINS
SUBROUTINE LTINV(KM,KMLOC,PSPVOR,PSPDIV,PSPSCALAR,FSPGL_PROC)

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS

USE PRLE1_MOD
USE PREPSNM_MOD
USE PRFI1_MOD
USE VDTUV_MOD
USE SPNSDE_MOD
USE LEINV_MOD
USE ASRE1_MOD
USE FSPGL_INT_MOD

#ifdef DOC

!**** *LTINV* - Inverse Legendre transform

!     Purpose.
!     --------
!        Tranform from Laplace space to Fourier space, compute U and V
!        and north/south derivatives of state variables.

!**   Interface.
!     ----------
!        *CALL* *LTINV(...)

!        Explicit arguments :
!        --------------------
!          KM        - zonal wavenumber
!          KMLOC     - local zonal wavenumber
!          PSPVOR    - spectral vorticity
!          PSPDIV    - spectral divergence
!          PSPSCALAR - spectral scalar variables

!        Implicit arguments :  The Laplace arrays of the model.
!        --------------------  The values of the Legendre polynomials
!                              The grid point arrays of the model
!     Method.
!     -------

!     Externals.
!     ----------
!         PRLE1   - prepares the Legendre polonymials
!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI1   - prepares the spectral fields
!         VDTUV   - compute u and v from vorticity and divergence
!         SPNSDE  - compute north-south derivatives
!         LEINV   - Inverse Legendre transform
!         ASRE1   - recombination of symmetric/antisymmetric part

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LTINV in IFS CY22R1
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE


INTEGER_M, INTENT(IN) :: KM
INTEGER_M, INTENT(IN) :: KMLOC
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

REAL_B  ,OPTIONAL, INTENT(IN) :: PSPVOR(:,:)
REAL_B  ,OPTIONAL, INTENT(IN) :: PSPDIV(:,:)
REAL_B  ,OPTIONAL, INTENT(IN) :: PSPSCALAR(:,:)

REAL_B :: ZIA(R%NLEI1,NLEI2)
REAL_B :: ZLEPO(R%NLEI3,R%NSMAX+2)
REAL_B :: ZEPSNM(0:R%NTMAX+2)
REAL_B :: ZSOA1(NLEI2,R%NLEI3),       ZAOA1(NLEI2,R%NLEI3)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IFC, ISTA
INTEGER_M :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISU,IDL,IDU

!     LOCAL LOGICAL SCALARS

!     LOCAL REAL SCALARS

!     ------------------------------------------------------------------

!*       1.    PREPARE LEGENDRE POLONOMIALS AND ZEPSNM.
!              ----------------------------------------

CALL PRLE1(KM,ZLEPO)

CALL PREPSNM(KM,KMLOC,ZEPSNM)

!     ------------------------------------------------------------------

!*       2.    PREPARE SPECTRAL FIELDS.
!              ------------------------

CALL PRFI1(KM,ZIA,PSPVOR,PSPDIV,PSPSCALAR)

!     ------------------------------------------------------------------

!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

IF (LUV) THEN
  IVORL = 1
  IVORU = 2*NF_UV
  IDIVL = 2*NF_UV+1
  IDIVU = 4*NF_UV
  IUL   = 4*NF_UV+1
  IUU   = 6*NF_UV
  IVL   = 6*NF_UV+1
  IVU   = 8*NF_UV
  CALL VDTUV(KM,NF_UV,ZEPSNM,ZIA(:,IVORL:IVORU),ZIA(:,IDIVL:IDIVU),&
             ZIA(:,IUL:IUU),ZIA(:,IVL:IVU))
ENDIF

IF (NF_SCDERS > 0) THEN
  ISL = 2*(4*NF_UV)+1
  ISU = ISL+2*NF_SCALARS-1
  IDL = 2*(4*NF_UV+NF_SCALARS)+1
  IDU = IDL+2*NF_SCDERS-1
  CALL SPNSDE(KM,ZEPSNM,ZIA(:,ISL:ISU),ZIA(:,IDL:IDU))
ENDIF

!     ------------------------------------------------------------------


!*       4.    INVERSE LEGENDRE TRANSFORM.
!              ---------------------------


ISTA = 1
IFC  = 2*NF_OUT_LT
IF(NF_UV > 0 .AND. .NOT. LVORGP) THEN
  ISTA = ISTA+2*NF_UV
  IFC  = IFC-2*NF_UV
ENDIF 
IF(NF_UV > 0 .AND. .NOT. LDIVGP) THEN
  ISTA = ISTA+2*NF_UV
  IFC  = IFC-2*NF_UV
ENDIF 

CALL LEINV(KM,IFC,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1(ISTA:,:),ZSOA1(ISTA:,:),ZLEPO)

!     ------------------------------------------------------------------

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
!              --------------------------------------------

CALL ASRE1(KM,KMLOC,ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!     6. OPTIONAL COMPUTATIONS IN FOURIER SPACE

IF(PRESENT(FSPGL_PROC)) THEN
  CALL FSPGL_INT(KM,KMLOC,FSPGL_PROC)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LTINV
END MODULE LTINV_MOD




