! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ELTINV_MOD
CONTAINS
SUBROUTINE ELTINV(KM,KMLOC,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B
USE TPMALD_DIM      ,ONLY : RALD
USE EPRFI1B_MOD     ,ONLY : EPRFI1B
USE EVDTUV_MOD      ,ONLY : EVDTUV
USE ESPNSDE_MOD     ,ONLY : ESPNSDE
USE ELEINV_MOD      ,ONLY : ELEINV
USE EASRE1B_MOD     ,ONLY : EASRE1B
USE FSPGL_INT_MOD   ,ONLY : FSPGL_INT
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

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
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        R. El Khatib 26-Aug-2021 Optimization for EASRE1B
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2
INTEGER(KIND=JPIM), INTENT(IN) :: KDIM1

REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPVOR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPDIV(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN) :: PSPMEANU(:)
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN) :: PSPMEANV(:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

REAL(KIND=JPRB) :: ZIA(RALD%NDGLSUR+R%NNOEXTZG,KLEI2)

INTEGER(KIND=JPIM) :: IFC, ISTA
INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISU,IDL,IDU
INTEGER(KIND=JPIM) :: IFIRST, ILAST,IDIM1,IDIM3,J3
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE



!     ------------------------------------------------------------------

!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTINV_MOD:ELTINV',0,ZHOOK_HANDLE)
IFIRST = 1
ILAST  = 4*KF_UV
ZIA=0.0_JPRB
IF (KF_UV > 0) THEN
  IVORL = 1
  IVORU = 2*KF_UV
  IDIVL = 2*KF_UV+1
  IDIVU = 4*KF_UV
  IUL   = 4*KF_UV+1
  IUU   = 6*KF_UV
  IVL   = 6*KF_UV+1
  IVU   = 8*KF_UV
  CALL EPRFI1B(KM,ZIA(:,IVORL:IVORU),PSPVOR,KF_UV,KFLDPTRUV)
  CALL EPRFI1B(KM,ZIA(:,IDIVL:IDIVU),PSPDIV,KF_UV,KFLDPTRUV)
  ILAST = ILAST+4*KF_UV
  CALL EVDTUV(KM,KF_UV,KFLDPTRUV,ZIA(:,IVORL:IVORU),ZIA(:,IDIVL:IDIVU),&
   & ZIA(:,IUL:IUU),ZIA(:,IVL:IVU),PSPMEANU,PSPMEANV)

ENDIF

IF(KF_SCALARS > 0)THEN
  IF(PRESENT(PSPSCALAR)) THEN
    IFIRST = ILAST+1
    ILAST  = IFIRST - 1 + 2*KF_SCALARS
    CALL EPRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
  ELSE
    IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
      IFIRST = ILAST+1
      ILAST  = IFIRST-1+2*NF_SC2
      CALL EPRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSC2(:,:),NF_SC2)
    ENDIF
    IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
      IDIM1=NF_SC3A
      IDIM3=UBOUND(PSPSC3A,3)
      DO J3=1,IDIM3
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*IDIM1
        CALL EPRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSC3A(:,:,J3),IDIM1)
      ENDDO
    ENDIF
    IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
      IDIM1=NF_SC3B
      IDIM3=UBOUND(PSPSC3B,3)
      DO J3=1,IDIM3
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*IDIM1
        CALL EPRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSC3B(:,:,J3),IDIM1)
      ENDDO
    ENDIF
  ENDIF
  IF(ILAST /= 8*KF_UV+2*KF_SCALARS) THEN
    WRITE(0,*) 'LTINV:KF_UV,KF_SCALARS,ILAST ',KF_UV,KF_SCALARS,ILAST
    CALL ABORT_TRANS('LTINV_MOD:ILAST /= 8*KF_UV+2*KF_SCALARS')
  ENDIF
ENDIF

IF (KF_SCDERS > 0) THEN
  ISL = 2*(4*KF_UV)+1
  ISU = ISL+2*KF_SCALARS-1
  IDL = 2*(4*KF_UV+KF_SCALARS)+1
  IDU = IDL+2*KF_SCDERS-1
  CALL ESPNSDE(KM,KF_SCALARS,ZIA(:,ISL:ISU),ZIA(:,IDL:IDU))
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

CALL ELEINV(KM,IFC,KF_OUT_LT,ZIA(:,ISTA:ISTA+IFC-1))

!     ------------------------------------------------------------------

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
!              --------------------------------------------

CALL EASRE1B(IFC,KM,KMLOC,ZIA(:,ISTA:ISTA+IFC-1))
!     ------------------------------------------------------------------

!     6. OPTIONAL COMPUTATIONS IN FOURIER SPACE

IF(PRESENT(FSPGL_PROC)) THEN
  CALL FSPGL_INT(KM,KMLOC,KF_UV,KF_SCALARS,KF_SCDERS,KF_OUT_LT,FSPGL_PROC,&
   & KFLDPTRUV,KFLDPTRSC)
ENDIF
IF (LHOOK) CALL DR_HOOK('ELTINV_MOD:ELTINV',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ELTINV
END MODULE ELTINV_MOD

