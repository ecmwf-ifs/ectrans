! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV_MOD
CONTAINS
SUBROUTINE LTINV(KM,KMLOC,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_FIELDS      ,ONLY : F
USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B, LATLON
USE TPM_FLT         ,ONLY : S
USE TPM_GEOMETRY    ,ONLY : G 

!USE PRLE1_MOD
USE PREPSNM_MOD     ,ONLY : PREPSNM
USE PRFI1B_MOD      ,ONLY : PRFI1B
USE VDTUV_MOD       ,ONLY : VDTUV
USE SPNSDE_MOD      ,ONLY : SPNSDE
USE LEINV_MOD       ,ONLY : LEINV
USE ASRE1B_MOD      ,ONLY : ASRE1B
USE FSPGL_INT_MOD   ,ONLY : FSPGL_INT
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE CDMAP_MOD       ,ONLY : CDMAP



!**** *LTINV* - Inverse Legendre transform
!
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
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

REAL(KIND=JPRB) :: ZACTHE
REAL(KIND=JPRB) :: ZIA(R%NLEI1,KLEI2)
REAL(KIND=JPRB) :: ZEPSNM(0:R%NTMAX+2)
!REAL(KIND=JPRB) :: ZSOA1(KDIM1,R%NLEI3),ZAOA1(KDIM1,R%NLEI3)
REAL(KIND=JPRB), ALLOCATABLE :: ZSOA1(:,:), ZAOA1(:,:), ZALN(:,:)

INTEGER(KIND=JPIM) :: IFC, ISTA, IIFC, IDGLU, JGL, JFLD
INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISLO,ISU,IDL,IDU, IGLS
INTEGER(KIND=JPIM) :: IFIRST, ILAST, IDIM1,IDIM3,J3
INTEGER(KIND=JPIM) :: INSDS, INSDE, IUVS, IUVE, IST

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!CHARACTER(LEN=10) :: CLHOOK

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!WRITE(CLHOOK,FMT='(A,I4.4)') 'LTINV_',KM
IF (LHOOK) CALL DR_HOOK('LTINV_MOD',0,ZHOOK_HANDLE)

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
           & ZIA(:,IUL:IUU),ZIA(:,IVL:IVU))
ENDIF

IF(KF_SCALARS > 0)THEN
  IF(PRESENT(PSPSCALAR)) THEN
    IFIRST = ILAST+1
    ILAST  = IFIRST - 1 + 2*KF_SCALARS
    CALL PRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
  ELSE
    IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
      IFIRST = ILAST+1
      ILAST  = IFIRST-1+2*NF_SC2
      CALL PRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSC2(:,:),NF_SC2)
    ENDIF
    IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
      IDIM1=NF_SC3A
      IDIM3=UBOUND(PSPSC3A,3)
      DO J3=1,IDIM3
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*IDIM1
        CALL PRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSC3A(:,:,J3),IDIM1)
      ENDDO
    ENDIF
    IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
      IDIM1=NF_SC3B
      IDIM3=UBOUND(PSPSC3B,3)
      DO J3=1,IDIM3
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*IDIM1
        CALL PRFI1B(KM,ZIA(:,IFIRST:ILAST),PSPSC3B(:,:,J3),IDIM1)
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

IIFC=IFC
IF(KM == 0)THEN
  IIFC=IFC/2
ENDIF

IF( LATLON.AND.S%LDLL ) THEN

  IDGLU = MIN(R%NDGNH,G%NDGLU(KM))

  IF( (S%LSHIFTLL .AND. KM < 2*IDGLU) .OR.&
   & (.NOT.S%LSHIFTLL .AND. KM < 2*(IDGLU-1)) ) THEN

    ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
    ISLO = S%FA(KMLOC)%ISLD

    ALLOCATE(ZAOA1(KDIM1,R%NLEI3))
    ALLOCATE(ZSOA1(KDIM1,R%NLEI3))
    CALL LEINV(KM,KMLOC,IFC,IIFC,KF_OUT_LT,ISL,IDGLU,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1,ZSOA1)

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
! before (non-linear) mapping !!!!

    ALLOCATE( ZALN(KDIM1,2*R%NDGNH) )
    DO JGL=ISL, R%NDGNH
      IGLS = 2*R%NDGNH+1-JGL
      DO JFLD=1,2*KF_OUT_LT
        ZALN(JFLD, JGL)  =  ZSOA1(JFLD,JGL)+ZAOA1(JFLD,JGL)
        ZALN(JFLD, IGLS) =  ZSOA1(JFLD,JGL)-ZAOA1(JFLD,JGL)
      ENDDO
    ENDDO
    
    IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
      IST = 1
      IF(LVORGP) THEN
        IST = IST+2*KF_UV
      ENDIF
      IF(LDIVGP) THEN
        IST = IST+2*KF_UV
      ENDIF
      IUVS = IST
      IUVE = IST+4*KF_UV-1
      IST = IST+4*KF_UV
      IST = IST+2*KF_SCALARS
      INSDS = IST
      INSDE = IST+2*KF_SCDERS-1
      IST = IST+2*KF_SCDERS
      
      IGLS = 2*R%NDGNH - ISL + 1
      IF( KF_UV > 0 ) THEN
        DO JGL=ISL, IGLS
          ZACTHE = REAL(F%RACTHE(JGL),JPRB)
          DO JFLD=IUVS,IUVE
            ZALN(JFLD, JGL)  =  ZALN(JFLD,JGL)*ZACTHE
          ENDDO
        ENDDO
      ENDIF
      IF( KF_SCDERS > 0 ) THEN
        DO JGL=ISL, IGLS
          ZACTHE = REAL(F%RACTHE(JGL),JPRB)
          DO JFLD=INSDS,INSDE
            ZALN(JFLD, JGL)  =  ZALN(JFLD,JGL)*ZACTHE
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    DEALLOCATE(ZAOA1)
    DEALLOCATE(ZSOA1)
    
    ! this routine maps to the output latitudes AND fills the FOUBUF
    CALL CDMAP(KM,KMLOC,ISL,ISLO,ZEPSNM(R%NTMAX+1),-1_JPIM,&
     & R%NDGNH,S%NDGNHD,2*KF_OUT_LT,ZALN,ZALN)
    DEALLOCATE(ZALN)
    
  ENDIF

ELSE
  IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
  ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

  ALLOCATE(ZAOA1(KDIM1,R%NLEI3))
  ALLOCATE(ZSOA1(KDIM1,R%NLEI3))
  CALL LEINV(KM,KMLOC,IFC,IIFC,KF_OUT_LT,ISL,IDGLU,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART/FILL FOUBUF
!              --------------------------------------------

  CALL ASRE1B(KF_OUT_LT,KM,KMLOC,ZAOA1,ZSOA1)
  DEALLOCATE(ZAOA1)
  DEALLOCATE(ZSOA1)
  
ENDIF

!     ------------------------------------------------------------------

!     6. OPTIONAL COMPUTATIONS IN FOURIER SPACE

IF(PRESENT(FSPGL_PROC)) THEN
  CALL FSPGL_INT(KM,KMLOC,KF_UV,KF_SCALARS,KF_SCDERS,KF_OUT_LT,FSPGL_PROC,&
   & KFLDPTRUV,KFLDPTRSC)
ENDIF

IF (LHOOK) CALL DR_HOOK('LTINV_MOD',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LTINV
END MODULE LTINV_MOD




