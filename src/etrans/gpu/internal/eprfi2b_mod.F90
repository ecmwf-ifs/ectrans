MODULE EPRFI2B_MOD
CONTAINS
SUBROUTINE EPRFI2B(KFIELD,PFFT,FOUBUF)

!**** *EPRFI2B* - Prepare input work arrays for direct transform

!     Purpose.
!     --------
!        To extract the Fourier fields for a specific zonal wavenumber
!        and put them in an order suitable for the direct Legendre
!        tranforms, i.e. split into symmetric and anti-symmetric part.

!**   Interface.
!     ----------
!     *CALL* *EPRFI2B(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM

!        Implicit arguments :  FOUBUF in TPM_TRANS
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
!        Original : 90-07-01
!        MPP Group: 95-10-01 Support for Distributed Memory version
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARKIND_ECTRANS, ONLY : JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_DISTR       ,ONLY : D
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
REAL(KIND=JPRBT)  , INTENT(OUT) :: PFFT(:,:,:)
REAL(KIND=JPRBT)  , INTENT(IN) :: FOUBUF(:)

INTEGER(KIND=JPIM) :: IM, JM
INTEGER(KIND=JPIM) :: ISTAN, JF, JGL
INTEGER(KIND=JPIM) :: IJR, IJI
REAL(KIND=JPRB) :: SCAL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI2B_MOD:EPRFI2B',0,ZHOOK_HANDLE)

!*       1.    EXTRACT SYM./ANTISYM. FIELDS FROM FOURIER ARRAY.
!              ------------------------------------------------

SCAL=1._JPRB/REAL(R_NDGL,JPRB)

!$acc data &
!$acc& present(PFFT) &
!$acc& present(FOUBUF) &
!$acc& copyin(R%NDGL,D%NPNTGTB1,D%NPROCL,D%NUMP,D%MYMS,SCAL)
  
!loop over wavenumber
!$acc parallel loop collapse(3) private(ISTAN,IM,IJR,IJI,JM)
DO JM = 1, D%NUMP
  DO JF =1,KFIELD
    DO JGL=1,R%NDGL
      IM = D%MYMS(JM)
      IJR = 2*(JF-1)+1
      IJI = IJR+1
      ISTAN = (D%NPNTGTB1(JM,JGL))*2*KFIELD
      PFFT(JGL,JM,IJR) = SCAL*FOUBUF(ISTAN+IJR)
      PFFT(JGL,JM,IJI) = SCAL*FOUBUF(ISTAN+IJI)
    ENDDO
  ENDDO
ENDDO

!$acc end data

IF (LHOOK) CALL DR_HOOK('EPRFI2B_MOD:EPRFI2B',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI2B
END MODULE EPRFI2B_MOD