MODULE LTINV_MOD
CONTAINS
SUBROUTINE LTINV(KM,KMLOC,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1,&
 & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS

USE PRLE1_MOD
USE PREPSNM_MOD
USE PRFI1B_MOD
USE VDTUV_MOD
USE SPNSDE_MOD
USE LEINV_MOD
USE ASRE1B_MOD
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

!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI1B  - prepares the spectral fields
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
INTEGER_M, INTENT(IN) :: KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1

REAL_B  ,OPTIONAL, INTENT(IN)  :: PSPVOR(:,:)
REAL_B  ,OPTIONAL, INTENT(IN)  :: PSPDIV(:,:)
REAL_B  ,OPTIONAL, INTENT(IN)  :: PSPSCALAR(:,:)
INTEGER_M,OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC 

REAL_B :: ZIA(R%NLEI1,KLEI2)
REAL_B :: ZEPSNM(0:R%NTMAX+2)
REAL_B :: ZSOA1(KDIM1,R%NLEI3),ZAOA1(KDIM1,R%NLEI3)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IFC, ISTA
INTEGER_M :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISU,IDL,IDU
INTEGER_M :: IDIV, IFIRST, ILAST, IVOR

!     LOCAL LOGICAL SCALARS

!     LOCAL REAL SCALARS

!     ------------------------------------------------------------------

!*       1.    PREPARE ZEPSNM.
!              ---------------

CALL PREPSNM(KM,KMLOC,ZEPSNM)

!     ------------------------------------------------------------------


!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

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
  CALL PRFI1B(KM,ZIA(:,IVORL:IVORU),PSPVOR,KF_UV,KFLDPTRUV)
  CALL PRFI1B(KM,ZIA(:,IDIVL:IDIVU),PSPDIV,KF_UV,KFLDPTRUV)
  ILAST = ILAST+4*KF_UV

  CALL VDTUV(KM,KF_UV,ZEPSNM,ZIA(:,IVORL:IVORU),ZIA(:,IDIVL:IDIVU),&
             ZIA(:,IUL:IUU),ZIA(:,IVL:IVU))
ENDIF

IF(KF_SCALARS > 0)THEN
  IFIRST = ILAST+1
  ILAST  = IFIRST - 1 + 2*KF_SCALARS 
  CALL PRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
ENDIF

IF (KF_SCDERS > 0) THEN
  ISL = 2*(4*KF_UV)+1
  ISU = ISL+2*KF_SCALARS-1
  IDL = 2*(4*KF_UV+KF_SCALARS)+1
  IDU = IDL+2*KF_SCDERS-1
  CALL SPNSDE(KM,KF_SCALARS,ZEPSNM,ZIA(:,ISL:ISU),ZIA(:,IDL:IDU))
ENDIF

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

CALL LEINV(KM,IFC,KF_OUT_LT,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
!              --------------------------------------------

CALL ASRE1B(KF_OUT_LT,KM,KMLOC,ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!     6. OPTIONAL COMPUTATIONS IN FOURIER SPACE

IF(PRESENT(FSPGL_PROC)) THEN
  CALL FSPGL_INT(KM,KMLOC,KF_UV,KF_SCALARS,KF_SCDERS,KF_OUT_LT,FSPGL_PROC,&
   & KFLDPTRUV,KFLDPTRSC)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LTINV
END MODULE LTINV_MOD




