MODULE EPRFI1AD_MOD
CONTAINS
SUBROUTINE EPRFI1AD(KM,KF_UV,KF_SCALARS,PIA,PSPVOR,PSPDIV,PSPSCALAR,&
 & KFLDPTRUV,KFLDPTRSC)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DISTR
!USE TPM_TRANS

USE EPRFI1BAD_MOD   ,ONLY : EPRFI1BAD

!**** *EPRFI1AD* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *EPRFI1AD(KM,PIA,PSPVOR,PSPDIV,PSPSCALAR

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
!        Original : 00-02-01 From PRFI1AD in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------
!
IMPLICIT NONE
!
!
INTEGER(KIND=JPIM),INTENT(IN) :: KM,KF_UV,KF_SCALARS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT)  :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT)  :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(INOUT)  :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,          INTENT(IN)     :: PIA(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)

INTEGER(KIND=JPIM) :: IDIV, IFIRST, ILAST, IVOR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              ------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI1AD_MOD:EPRFI1AD',0,ZHOOK_HANDLE)
IFIRST = 1
ILAST  = 4*KF_UV

!*       1.1      VORTICITY AND DIVERGENCE.

IF(KF_UV > 0)THEN
  IVOR = 1
  IDIV = 2*KF_UV+1
  CALL EPRFI1BAD(KM,PIA(:,IVOR:IDIV-1),PSPVOR,KF_UV,KFLDPTRUV)
  CALL EPRFI1BAD(KM,PIA(:,IDIV:ILAST) ,PSPDIV,KF_UV,KFLDPTRUV)
  ILAST = ILAST+4*KF_UV
ENDIF

!*       1.2    SCALAR VARIABLES.

IF(KF_SCALARS > 0)THEN
  IFIRST = ILAST+1
  ILAST  = IFIRST - 1 + 2*KF_SCALARS
  CALL EPRFI1BAD(KM,PIA(:,IFIRST:ILAST),PSPSCALAR(:,:),KF_SCALARS,KFLDPTRSC)
ENDIF
IF (LHOOK) CALL DR_HOOK('EPRFI1AD_MOD:EPRFI1AD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI1AD
END MODULE EPRFI1AD_MOD
