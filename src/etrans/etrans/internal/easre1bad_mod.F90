MODULE EASRE1BAD_MOD
CONTAINS
SUBROUTINE EASRE1BAD(KFC,KM,KMLOC,PIA)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPMALD_DIM      ,ONLY : RALD
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_DISTR       ,ONLY : D

!**** *EASRE1BAD* - Recombine antisymmetric and symmetric parts - adjoint

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!        *CALL* *EASRE1BAD(..)

!        Explicit arguments :
!        -------------------   KFC - number of fields (input-c)
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
!        Original : 00-02-01 From ASRE1BAD in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib 26-Aug-2021 Optimizations
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFC
INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC

REAL(KIND=JPRB), INTENT(OUT)    :: PIA(RALD%NDGLSUR+R%NNOEXTZG,KFC)

INTEGER(KIND=JPIM) ::   JFLD, JGL ,IPROC
INTEGER(KIND=JPIM) :: IISTAN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    RECOMBINATION  OF SYMMETRIC AND ANTSYMMETRIC PARTS.
!              ---------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EASRE1BAD_MOD:EASRE1BAD',0,ZHOOK_HANDLE)
#ifdef __INTEL_COMPILER
!$OMP SIMD PRIVATE(JGL)
DO JFLD  =1,KFC
  DO JGL=1,R%NDGL
    PIA(JGL,JFLD)=FOUBUF_IN((D%NSTAGT0B(D%NPROCL(JGL))+D%NPNTGTB1(KMLOC,JGL))*KFC+JFLD)
  ENDDO
ENDDO
#else
DO JGL=1,R%NDGL
  IPROC=D%NPROCL(JGL)
  DO JFLD  =1,KFC
    IISTAN=(D%NSTAGT0B(IPROC) + D%NPNTGTB1(KMLOC,JGL))*KFC
    PIA(JGL,JFLD)=FOUBUF_IN(IISTAN+JFLD)
  ENDDO
ENDDO
#endif
IF (LHOOK) CALL DR_HOOK('EASRE1BAD_MOD:EASRE1BAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EASRE1BAD
END MODULE EASRE1BAD_MOD
