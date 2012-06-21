MODULE LTINVAD_MOD
CONTAINS
SUBROUTINE LTINVAD(KM,KMLOC,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1,&
 & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS

USE PRLE1AD_MOD
USE PREPSNM_MOD
USE PRFI1BAD_MOD
USE VDTUVAD_MOD
USE SPNSDEAD_MOD
USE LEINVAD_MOD
USE ASRE1BAD_MOD
USE FSPGL_INT_MOD

#ifdef DOC

!**** *LTINVAD* - Inverse Legendre transform

!     Purpose.
!     --------
!        Tranform from Laplace space to Fourier space, compute U and V
!        and north/south derivatives of state variables.

!**   Interface.
!     ----------
!        *CALL* *LTINVAD(...)

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
!         PRLE1AD - prepares the Legendre polonymials
!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI1AD - prepares the spectral fields
!         VDTUVAD - compute u and v from vorticity and divergence
!         SPNSDEAD- compute north-south derivatives
!         LEINVAD - Inverse Legendre transform
!         ASRE1AD - recombination of symmetric/antisymmetric part

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LTINVAD in IFS CY22R1
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE


INTEGER_M, INTENT(IN) :: KM
INTEGER_M, INTENT(IN) :: KMLOC
INTEGER_M, INTENT(IN) :: KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1

REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
INTEGER_M,OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

REAL_B :: ZIA(R%NLEI1,KLEI2)
REAL_B :: ZEPSNM(0:R%NTMAX+2)
REAL_B :: ZSOA1(KDIM1,R%NLEI3),ZAOA1(KDIM1,R%NLEI3)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IFC, ISTA
INTEGER_M :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISU,IDL,IDU
INTEGER_M :: ILAST,IFIRST

!     LOCAL LOGICAL SCALARS

!     LOCAL REAL SCALARS

!     ------------------------------------------------------------------

!*       1.    PREPARE AND ZEPSNM.
!              -------------------

CALL PREPSNM(KM,KMLOC,ZEPSNM)

!     ------------------------------------------------------------------
!     6. OPTIONAL COMPUTATIONS IN FOURIER SPACE

!IF(PRESENT(FSPGL_PROC)) THEN
!  CALL FSPGL_INT(KM,KMLOC,FSPGL_PROC)
!ENDIF

!     ------------------------------------------------------------------

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
!              --------------------------------------------

CALL ASRE1BAD(KF_OUT_LT,KM,KMLOC,ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!*       4.    INVERSE LEGENDRE TRANSFORM.
!              ---------------------------


ISTA = 1
IFC  = 2*KF_OUT_LT
IF(KF_UV > 0 .AND. .NOT. LVORGP) THEN
  ISTA = ISTA+2*KF_UV
ENDIF 
IF(KF_UV > 0 .AND. .NOT. LDIVGP) THEN
  ISTA = ISTA+2*KF_UV
ENDIF 

ZIA(:,ISTA:ISTA+IFC-1) = _ZERO_

CALL LEINVAD(KM,IFC,KF_OUT_LT,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

ZIA(:,1:ISTA-1) = _ZERO_

IFIRST = 1
ILAST  = 4*KF_UV
IF (KF_UV > 0) THEN
  IVORL = 1
  IVORU = 2*KF_UV
  IDIVL = 2*KF_UV+1
  IDIVU = 4*KF_UV
  IUL   = 4*KF_UV+1
  IUU   = 6*KF_UV
  IVL   = 6*KF_UV+1
  IVU   = 8*KF_UV
  CALL VDTUVAD(KM,KF_UV,ZEPSNM,ZIA(:,IVORL:IVORU),ZIA(:,IDIVL:IDIVU),&
             ZIA(:,IUL:IUU),ZIA(:,IVL:IVU))
  CALL PRFI1BAD(KM,ZIA(:,IVORL:IVORU),PSPVOR,KF_UV,KFLDPTRUV)
  CALL PRFI1BAD(KM,ZIA(:,IDIVL:IDIVU),PSPDIV,KF_UV,KFLDPTRUV)
  ILAST = ILAST+4*KF_UV
ENDIF

IF (KF_SCDERS > 0) THEN
  ISL = 2*(4*KF_UV)+1
  ISU = ISL+2*KF_SCALARS-1
  IDL = 2*(4*KF_UV+KF_SCALARS)+1
  IDU = IDL+2*KF_SCDERS-1
  CALL SPNSDEAD(KM,KF_SCALARS,ZEPSNM,ZIA(:,ISL:ISU),ZIA(:,IDL:IDU))
ENDIF

IF(KF_SCALARS > 0)THEN
  IFIRST = ILAST+1
  ILAST  = IFIRST - 1 + 2*KF_SCALARS 
  CALL PRFI1BAD(KM,ZIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
ENDIF



!     ------------------------------------------------------------------


END SUBROUTINE LTINVAD
END MODULE LTINVAD_MOD




