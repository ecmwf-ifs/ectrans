MODULE EUPDSP_MOD
CONTAINS
SUBROUTINE EUPDSP(KM,KF_UV,KF_SCALARS,PFFT,PVODI, &
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC)

!**** *EUPDSP* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update the spectral arrays for a fixed zonal wave-number
!        from values in POA1 and POA2.

!**   Interface.
!     ----------
!        CALL EUPDSP(...)

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
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_TRANS       ,ONLY : NF_SC2, NF_SC3A, NF_SC3B
!USE TPM_DISTR

USE EUPDSPB_MOD     ,ONLY : EUPDSPB
!

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN)  :: KM,KF_UV,KF_SCALARS
REAL(KIND=JPRB) , INTENT(IN)  :: PFFT(:,:)
REAL(KIND=JPRB) , INTENT(IN)  :: PVODI(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)

INTEGER(KIND=JPIM) :: IVORS, IVORE, IDIVS, IDIVE, IST ,IEND,IDIM1,IDIM3,J3
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    UPDATE FIELDS
!              -------------

!*       1.1      VORTICITY AND DIVERGENCE.

IF (LHOOK) CALL DR_HOOK('EUPDSP_MOD:EUPDSP',0,ZHOOK_HANDLE)
IST = 1
IF (KF_UV > 0) THEN
  IST = IST+4*KF_UV
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV
  CALL EUPDSPB(KM,KF_UV,PVODI(:,IVORS:IVORE),PSPVOR,KFLDPTRUV)
  CALL EUPDSPB(KM,KF_UV,PVODI(:,IDIVS:IDIVE),PSPDIV,KFLDPTRUV)
ENDIF

!*       1.2   SCALARS

IF (KF_SCALARS > 0) THEN
  IF(PRESENT(PSPSCALAR)) THEN
    IEND = IST+2*KF_SCALARS-1
    CALL EUPDSPB(KM,KF_SCALARS,PFFT(:,IST:IEND),PSPSCALAR,KFLDPTRSC)
  ELSE
    IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
      IDIM1 = NF_SC2
      IEND  = IST+2*IDIM1-1
      CALL EUPDSPB(KM,IDIM1,PFFT(:,IST:IEND),PSPSC2)
      IST=IST+2*IDIM1
    ENDIF
    IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
      IDIM1=NF_SC3A
      IDIM3=UBOUND(PSPSC3A,3)
      DO J3=1,IDIM3
        IEND = IST+2*IDIM1-1
        CALL EUPDSPB(KM,IDIM1,PFFT(:,IST:IEND),PSPSC3A(:,:,J3))
        IST=IST+2*IDIM1
      ENDDO
    ENDIF
    IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
      IDIM1=NF_SC3B
      IDIM3=UBOUND(PSPSC3B,3)
      DO J3=1,IDIM3
        IEND = IST+2*IDIM1-1
        CALL EUPDSPB(KM,IDIM1,PFFT(:,IST:IEND),PSPSC3B(:,:,J3))
        IST=IST+2*IDIM1
      ENDDO
    ENDIF
  ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('EUPDSP_MOD:EUPDSP',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EUPDSP
END MODULE EUPDSP_MOD
