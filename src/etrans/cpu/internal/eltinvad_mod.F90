! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ELTINVAD_MOD
CONTAINS
SUBROUTINE ELTINVAD(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,&
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC,PSPMEANU,PSPMEANV)

!**** *ELTINVAD* - Control routine for inverse Legandre transform - adj.

!     Purpose.
!     --------
!     Control routine for the inverse LEGENDRE transform

!**   Interface.
!     ----------
!     CALL ELTINVAD(...)
!     KF_OUT_LT    - number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KFLDPTRUV(:) - field pointer array for vor./div.
!     KFLDPTRSC(:) - field pointer array for PSPSCALAR
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition  

!     Method.
!     -------

!     Externals.  
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LTINVAD in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        01-Dec-2004   A. Deckmyn    add KMLOC to EVDTUVAD call
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement +
!        thread-safety
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN
USE TPM_DIM         ,ONLY : R
USE TPMALD_DIM      ,ONLY : RALD
USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B
USE TPM_DISTR

USE EASRE1BAD_MOD   ,ONLY : EASRE1BAD
USE ELEINVAD_MOD    ,ONLY : ELEINVAD
USE EPRFI1BAD_MOD   ,ONLY : EPRFI1BAD
USE ESPNSDEAD_MOD   ,ONLY : ESPNSDEAD
USE EVDTUVAD_MOD    ,ONLY : EVDTUVAD
USE EVDTUVAD_COMM_MOD
USE EXTPER_MOD      ,ONLY : EXTPER

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2

REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPMEANU(:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(INOUT) :: PSPMEANV(:)

INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)

EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

REAL(KIND=JPRB) :: ZIA(RALD%NDGLSUR+R%NNOEXTZG,KLEI2,D%NUMP)
REAL(KIND=JPRB) :: ZIA2(KLEI2,RALD%NDGLSUR+R%NNOEXTZG)

INTEGER(KIND=JPIM) :: IFC, ISTA, IINDEX(2*KF_OUT_LT), JF, JDIM, IM, JM
INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISU,IDL,IDU
INTEGER(KIND=JPIM) :: ILAST,IFIRST,IDIM1,IDIM3,J3

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTINVAD_MOD:ELTINVAD',0,ZHOOK_HANDLE)

IF (KF_UV > 0) THEN
  IVORL = 1
  IVORU = 2*KF_UV
  IDIVL = 2*KF_UV+1
  IDIVU = 4*KF_UV
  IUL   = 4*KF_UV+1
  IUU   = 6*KF_UV
  IVL   = 6*KF_UV+1
  IVU   = 8*KF_UV
ENDIF
ISTA = 1
IFC  = 2*KF_OUT_LT
IF(KF_UV > 0 .AND. .NOT. LVORGP) THEN
  ISTA = ISTA+2*KF_UV
ENDIF
IF(KF_UV > 0 .AND. .NOT. LDIVGP) THEN
  ISTA = ISTA+2*KF_UV
ENDIF
IF (KF_SCDERS > 0) THEN
  ISL = 2*(4*KF_UV)+1
  ISU = ISL+2*KF_SCALARS-1
  IDL = 2*(4*KF_UV+KF_SCALARS)+1
  IDU = IDL+2*KF_SCDERS-1
ENDIF

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM,JF,JDIM,IINDEX,ZIA2)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)

!        7.    OPTIONAL COMPUTATIONS IN FOURIER SPACE
!              --------------------------------------

!commented  IF(PRESENT(FSPGL_PROC)) THEN
!commented    CALL FSPGL_INT(IM,JM,FSPGL_PROC)
!commented  ENDIF


!*       6.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
!              --------------------------------------------

  ZIA(:,:,JM)=0.0_JPRB
  CALL EASRE1BAD(IFC,IM,JM,ZIA(:,ISTA:ISTA+IFC-1,JM))


!*       5.    PERIODICIZATION IN Y DIRECTION
!              ------------------------------

  IF(R%NNOEXTZG>0) THEN
    DO JF = 1,IFC
      DO JDIM = 1,R%NDGL
        ZIA2(JF,JDIM)=ZIA(JDIM,JF,JM)
      ENDDO
    ENDDO
    IINDEX(1)=0
    CALL EXTPER(ZIA2(:,:),R%NDGL+R%NNOEXTZG,1,R%NDGL,IFC,1,IINDEX,0)
    DO JF = 1,IFC
      DO JDIM = 1,R%NDGL+R%NNOEXTZG
        ZIA(JDIM,JF,JM) = ZIA2(JF,JDIM)
      ENDDO
    ENDDO
  ENDIF

!*       4.    INVERSE LEGENDRE TRANSFORM.
!              ---------------------------

  CALL ELEINVAD(IM,IFC,KF_OUT_LT,ZIA(:,ISTA:ISTA+IFC-1,JM))


!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

  ZIA(:,1:ISTA-1,JM) = 0.0_JPRB

  IF (KF_UV > 0) THEN
    CALL EVDTUVAD(IM,JM,KF_UV,KFLDPTRUV,ZIA(:,IVORL:IVORU,JM),ZIA(:,IDIVL:IDIVU,JM),&
     & ZIA(:,IUL:IUU,JM),ZIA(:,IVL:IVU,JM),PSPMEANU,PSPMEANV)
  ENDIF


ENDDO
!$OMP END PARALLEL DO

!*       2.    COMMUNICATION OF MEAN WIND
!              --------------------------

IF (KF_UV > 0) THEN
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL EVDTUVAD_COMM(IM,JM,KF_UV,KFLDPTRUV,PSPMEANU,PSPMEANV)
  ENDDO
ENDIF

!*       2.    PREPARE SPECTRAL FIELDS
!              -----------------------

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM,IFIRST,ILAST,IDIM1,IDIM3)
DO JM=1,D%NUMP
  IM = D%MYMS(JM)

  IFIRST = 1
  ILAST  = 4*KF_UV
  IF (KF_UV > 0) THEN
    CALL EPRFI1BAD(IM,ZIA(:,IVORL:IVORU,JM),PSPVOR,KF_UV,KFLDPTRUV)
    CALL EPRFI1BAD(IM,ZIA(:,IDIVL:IDIVU,JM),PSPDIV,KF_UV,KFLDPTRUV)
    ILAST = ILAST+4*KF_UV
  ENDIF

  IF (KF_SCDERS > 0) THEN
    CALL ESPNSDEAD(IM,KF_SCALARS,ZIA(:,ISL:ISU,JM),ZIA(:,IDL:IDU,JM))
  ENDIF

  IF(KF_SCALARS > 0)THEN
    IF(PRESENT(PSPSCALAR)) THEN
      IFIRST = ILAST+1
      ILAST  = IFIRST - 1 + 2*KF_SCALARS 
      CALL EPRFI1BAD(IM,ZIA(:,IFIRST:ILAST,JM),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
    ELSE
      IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*NF_SC2
        CALL EPRFI1BAD(IM,ZIA(:,IFIRST:ILAST,JM),PSPSC2(:,:),NF_SC2)
      ENDIF
      IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
        IDIM1=NF_SC3A
        IDIM3=UBOUND(PSPSC3A,3)
        DO J3=1,IDIM3
          IFIRST = ILAST+1
          ILAST  = IFIRST-1+2*IDIM1
          CALL EPRFI1BAD(IM,ZIA(:,IFIRST:ILAST,JM),PSPSC3A(:,:,J3),IDIM1)
        ENDDO
      ENDIF
      IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
        IDIM1=NF_SC3B
        IDIM3=UBOUND(PSPSC3B,3)
        DO J3=1,IDIM3
          IFIRST = ILAST+1
          ILAST  = IFIRST-1+2*IDIM1
          CALL EPRFI1BAD(IM,ZIA(:,IFIRST:ILAST,JM),PSPSC3B(:,:,J3),IDIM1)
        ENDDO
      ENDIF
    ENDIF
  ENDIF

ENDDO
!$OMP END PARALLEL DO

IF (LHOOK) CALL DR_HOOK('ELTINVAD_MOD:ELTINVAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ELTINVAD
END MODULE ELTINVAD_MOD
