MODULE EPRFI1B_MOD
CONTAINS
SUBROUTINE EPRFI1B(PFFT,PSPEC,KFIELDS,KFLDPTR)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DIM
USE TPM_DISTR
USE TPMALD_DISTR    ,ONLY : DALD, DALD_NESM0, DALD_NCPL2M
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

REAL(KIND=JPRB)   ,INTENT(OUT)  :: PFFT(:,:,:)
REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF,IFLD
INTEGER(KIND=JPIM) :: IM, JM, MAX_NCPL2M
INTEGER(KIND=JPIM) :: JFLDPTR(KFIELDS)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
character(len=64) :: frmt

!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              --------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',0,ZHOOK_HANDLE)

IF (PRESENT(KFLDPTR)) THEN
  JFLDPTR=KFLDPTR
ELSE
  DO JFLD=1,KFIELDS
    JFLDPTR(JFLD)=JFLD
  ENDDO
ENDIF

!$acc data present (PFFT, PSPEC)

!!$acc kernels default(none)
!PFFT = 0._JPRB
!!$acc end kernels

MAX_NCPL2M = MAXVAL (DALD_NCPL2M)

!$ACC parallel loop collapse(3) &
!$ACC& present(D_MYMS,DALD_NCPL2M,DALD_NESM0,D_NUMP) &
!$ACC& copyin(KFIELDS,MAX_NCPL2M,JFLDPTR) &
!$ACC& private(IR,II,IM,ILCM,IOFF,INM,JFLD) default(none)
DO JM = 1, D_NUMP
  DO JFLD=1,KFIELDS
    DO J=1,MAX_NCPL2M,2
      IR = 2*JFLD-1
      II = IR+1
      IM   = D_MYMS(JM)
      ILCM = DALD_NCPL2M(IM)
      IF (J .LE. ILCM) THEN
        IOFF = DALD_NESM0(IM)
        INM = IOFF+(J-1)*2
        PFFT(J  ,JM,IR) = PSPEC(JFLDPTR(JFLD),INM  )
        PFFT(J+1,JM,IR) = PSPEC(JFLDPTR(JFLD),INM+1)
        PFFT(J  ,JM,II) = PSPEC(JFLDPTR(JFLD),INM+2)
        PFFT(J+1,JM,II) = PSPEC(JFLDPTR(JFLD),INM+3)
      ELSE
        PFFT(J  ,JM,IR) = 0._JPRB
        PFFT(J+1,JM,IR) = 0._JPRB
        PFFT(J  ,JM,II) = 0._JPRB
        PFFT(J+1,JM,II) = 0._JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$acc end data



IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI1B
END MODULE EPRFI1B_MOD