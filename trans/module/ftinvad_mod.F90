MODULE FTINVAD_MOD
CONTAINS
SUBROUTINE FTINVAD(PREEL,KFIELDS,KGL)


!**** *FTINVAD - Inverse Fourier transform - adjoint

!     Purpose. Routine for Fourier to Grid-point transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINVAD(..)

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
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR
USE TPM_TRANS
USE TPM_GEOMETRY
USE TPM_FFT
USE TPM_DIM

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(OUT)  :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN,IOFF

!     ------------------------------------------------------------------

IJUMP = 1
IGLG  = D%NPTRLS(MYSETW)+KGL-1
ILOEN = G%NLOEN(IGLG)+R%NNOEXTZL
IST   = 2*(G%NMEN(IGLG)+1)+1
ILEN  = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)

  ! Change of metric (not in forward routine)
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF+JJ) = PREEL(JF,IOFF+JJ)*ILOEN
  ENDDO
ENDDO
  
CALL FFT992(PREEL(1,IOFF+1),T%TRIGS(1,KGL),&
 &T%NFAX(1,KGL),KFIELDS,IJUMP,ILOEN,KFIELDS,-1)

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+IOFF+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTINVAD
END MODULE FTINVAD_MOD
