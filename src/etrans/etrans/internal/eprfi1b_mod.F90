MODULE EPRFI1B_MOD
CONTAINS
SUBROUTINE EPRFI1B(KM,PIA,PSPEC,KFIELDS,KFLDPTR)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DIM
!USE TPM_DISTR
USE TPMALD_DISTR    ,ONLY : DALD
!
!**** *PRFI1* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1B(...)*

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              PSPEC  - spectral array
!                              KFIELDS  - number of fields

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
!        Original : 00-02-01 From PRFI1B in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KM,KFIELDS
REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)  :: PIA(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF,IFLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              --------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',0,ZHOOK_HANDLE)
ILCM = DALD%NCPL2M(KM)
IOFF = DALD%NESM0(KM)

IF(PRESENT(KFLDPTR)) THEN
  DO JFLD=1,KFIELDS
    IR = 2*(JFLD-1)+1
    II = IR+1
    IFLD = KFLDPTR(JFLD)
    DO J=1,ILCM,2
      INM = IOFF+(J-1)*2
      PIA(J  ,IR) = PSPEC(IFLD,INM  )
      PIA(J+1,IR) = PSPEC(IFLD,INM+1)
      PIA(J  ,II) = PSPEC(IFLD,INM+2)
      PIA(J+1,II) = PSPEC(IFLD,INM+3)
    ENDDO
  ENDDO

ELSE
  DO J=1,ILCM,2
    INM = IOFF+(J-1)*2
    !DIR$ IVDEP
    !OCL NOVREC
    !cdir unroll=4
    DO JFLD=1,KFIELDS
      IR = 2*(JFLD-1)+1
      II = IR+1
      PIA(J  ,IR) = PSPEC(JFLD,INM  )
      PIA(J+1,IR) = PSPEC(JFLD,INM+1)
      PIA(J  ,II) = PSPEC(JFLD,INM+2)
      PIA(J+1,II) = PSPEC(JFLD,INM+3)
    ENDDO
  ENDDO

ENDIF
IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI1B
END MODULE EPRFI1B_MOD
