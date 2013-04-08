MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(PREEL,KFIELDS,KGL)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

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
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR       ,ONLY : D, MYSETW
!USE TPM_TRANS
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFT         ,ONLY : T
USE TPM_DIM         ,ONLY : R
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1

!     ------------------------------------------------------------------

IJUMP = 1
IGLG = D%NPTRLS(MYSETW)+KGL-1
IST  = 2*(G%NMEN(IGLG)+1)+1
ILEN = G%NLOEN(IGLG)+R%NNOEXTZL+3-IST

IF (G%NLOEN(IGLG)>1) THEN 
  CALL FFT992(PREEL(1,D%NSTAGTF(KGL)+1),T%TRIGS(1,KGL),&
   &T%NFAX(1,KGL),KFIELDS,IJUMP,G%NLOEN(IGLG)+R%NNOEXTZL,KFIELDS,-1)
ENDIF
IST1=1
IF (G%NLOEN(IGLG)==1) IST1=0
DO JJ=IST1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+D%NSTAGTF(KGL)+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
