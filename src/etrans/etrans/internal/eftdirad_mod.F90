MODULE EFTDIRAD_MOD
CONTAINS
SUBROUTINE EFTDIRAD(PREEL,KFIELDS,KGL)

!**** *EFTDIRAD - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform - adjoint
!     --------

!**   Interface.
!     ----------
!        CALL EFTDIRAD(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib 01-Sep-2015 support for FFTW transforms
!        R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, MYSETW
!USE TPM_TRANS
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FFT         ,ONLY : T, TB
USE BLUESTEIN_MOD   ,ONLY : BLUESTEIN_FFT
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_DIM         ,ONLY : R
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS,KGL
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,ILOEN
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN,ITYPE
REAL(KIND=JPRB) :: ZNORM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EFTDIRAD_MOD:EFTDIRAD',0,ZHOOK_HANDLE)

ITYPE = 1
IJUMP = 1
IGLG = D%NPTRLS(MYSETW)+KGL-1
ILOEN = G%NLOEN(IGLG)
IST  = 2*(G%NMEN(IGLG)+1)+1
ILEN = ILOEN+3-IST
IOFF  = D%NSTAGTF(KGL)+1

DO JJ=1,ILEN
  DO JF=1,KFIELDS
    PREEL(JF,IST+IOFF-1+JJ-1) = 0.0_JPRB
  ENDDO
ENDDO
DO JJ=1,1
  DO JF=1,KFIELDS
    PREEL(JF,IOFF-1+JJ) = 2.0_JPRB * PREEL(JF,IOFF-1+JJ)
  ENDDO
ENDDO

#ifdef WITH_FFTW
IF( .NOT. TW%LFFTW )THEN
#endif

  IF( T%LUSEFFT992(KGL) )THEN

    CALL FFT992(PREEL(1,IOFF),T%TRIGS(1,KGL),&
     &T%NFAX(1,KGL),KFIELDS,IJUMP,ILOEN,KFIELDS,ITYPE)

  ELSE
  
    CALL BLUESTEIN_FFT(TB,ILOEN,ITYPE,KFIELDS,PREEL(1:KFIELDS,IOFF:IOFF+ICLEN-1))

  ENDIF

#ifdef WITH_FFTW
ELSE

  IRLEN=G%NLOEN(IGLG)+R%NNOEXTZL
  ICLEN=(IRLEN/2+1)*2
  CALL EXEC_FFTW(ITYPE,IRLEN,ICLEN,IOFF,KFIELDS,TW%LALL_FFTW,PREEL)

ENDIF
#endif


  ! Change of metric (not in forward routine)
ZNORM=1.0_JPRB/(2.0_JPRB*REAL(ILOEN,JPRB))
DO JJ=1,ILOEN
  DO JF=1,KFIELDS
    PREEL(JF,IOFF-1+JJ) = ZNORM * PREEL(JF,IOFF-1+JJ)
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('EFTDIRAD_MOD:EFTDIRAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EFTDIRAD
END MODULE EFTDIRAD_MOD
