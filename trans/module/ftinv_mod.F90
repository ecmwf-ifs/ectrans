MODULE FTINV_MOD
CONTAINS
SUBROUTINE FTINV(PREEL,KFIELDS,KGL)

!**** *FTINV - Inverse Fourier transform

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------  
!                 

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 : 2D model (NLOEN=1)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR
USE TPM_TRANS
USE TPM_GEOMETRY
USE TPM_FFT

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(OUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1

!     ------------------------------------------------------------------

IJUMP = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = G%NLOEN(IGLG)+3-IST
IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0

DO JJ=IST1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

IF (G%NLOEN(IGLG)>1) THEN
  CALL FFT992(PREEL(1,D%NSTAGTF(KGL)+1),T%TRIGS(1,KGL),&
   &T%NFAX(1,KGL),KFIELDS,IJUMP,G%NLOEN(IGLG),KFIELDS,1)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FTINV
END MODULE FTINV_MOD
