MODULE EASRE1B_MOD
CONTAINS
SUBROUTINE EASRE1B(KFIELD,PFFT,FOUBUF_IN)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D

!**** *ASRE1B* - Recombine antisymmetric and symmetric parts

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!        *CALL* *ASRE1B(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields (input-c)
!                              KM - zonal wavenumber(input-c)
!                              KMLOC - local version of KM (input-c)
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (input)
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM (input)

!        Implicit arguments :  FOUBUF_IN - output buffer (output)
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
!        Original : 00-02-01 From ASRE1B in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRB), INTENT(IN)    :: PFFT(:,:,:)
REAL(KIND=JPRB), INTENT(OUT) :: FOUBUF_IN(:)

INTEGER(KIND=JPIM) :: JFLD, JGL ,IPROC
INTEGER(KIND=JPIM) :: IISTAN
INTEGER(KIND=JPIM) :: JM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    RECOMBINATION  OF SYMMETRIC AND ANTSYMMETRIC PARTS.
!              ---------------------------------------------------


IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',0,ZHOOK_HANDLE)

!$acc parallel loop collapse(3) private (JM, JGL, JFLD, IPROC, IISTAN) &
!$acc& present (PFFT, D%NSTAGT0B, D%NPNTGTB1, D%NPROCL, D%NUMP, R%NDGL, FOUBUF_IN) &
!$acc& copyin(KFIELD) default(none)
DO JM = 1, D%NUMP  !100
  DO JGL=1,R%NDGL  !400
    DO JFLD  =1,2*KFIELD !500
      IPROC=D%NPROCL(JGL)
      IISTAN=(D%NPNTGTB1(JM,JGL))*2*KFIELD
      FOUBUF_IN(IISTAN+JFLD)=PFFT(JGL,JM,JFLD)
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

IF (LHOOK) CALL DR_HOOK('EASRE1B_MOD:EASRE1B',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE EASRE1B
END MODULE EASRE1B_MOD