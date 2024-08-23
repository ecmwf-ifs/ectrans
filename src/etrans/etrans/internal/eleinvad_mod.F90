MODULE ELEINVAD_MOD
CONTAINS
SUBROUTINE ELEINVAD(KM,KFC,KF_OUT_LT,PIA)

!**** *ELEINVAD* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL ELEINVAD(...)

!        Explicit arguments :  KM - zonal wavenumber (input-c)
!        --------------------  KFC - number of fields to tranform (input-c)
!                              PIA - spectral fields
!                              for zonal wavenumber KM (input)
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PLEPO - Legendre polonomials for zonal
!                              wavenumber KM (input-c)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - calls SGEMVX (matrix multiply)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LEINVAD in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        R. El Khatib 01-Sep-2015 support for FFTW transforms
!        R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
!USE TPM_GEOMETRY
!USE TPM_TRANS
USE TPMALD_FFT      ,ONLY : TALD
#ifdef WITH_FFTW
USE TPM_FFTW     ,ONLY : TW, EXEC_EFFTW
#endif
USE TPMALD_DIM      ,ONLY : RALD
USE TPMALD_FFT      ,ONLY : TALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KM
INTEGER(KIND=JPIM), INTENT(IN)    :: KFC
INTEGER(KIND=JPIM), INTENT(IN)    :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(OUT)   :: PIA(:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN, IOFF, ITYPE
INTEGER(KIND=JPIM) :: JJ, JF
REAL(KIND=JPRB) :: ZNORM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELEINVAD_MOD:ELEINVAD',0,ZHOOK_HANDLE)

IF (KFC>0) THEN
  ITYPE=-1
  IRLEN=R%NDGL+R%NNOEXTZG
  ICLEN=RALD%NDGLSUR+R%NNOEXTZG
  IF( TALD%LFFT992 )THEN
    CALL FFT992(PIA,TALD%TRIGSE,TALD%NFAXE,1,ICLEN,IRLEN,KFC,ITYPE)
#ifdef WITH_FFTW
  ELSEIF( TW%LFFTW )THEN
    IOFF=1
    CALL EXEC_EFFTW(ITYPE,IRLEN,ICLEN,IOFF,KFC,TW%LALL_FFTW,PIA)
#endif
  ELSE
    CALL ABORT_TRANS('ELEDIR_MOD:ELEINVAD: NO FFT PACKAGE SELECTED')
  ENDIF
  ZNORM=2.0_JPRB*REAL(R%NDGL+R%NNOEXTZG,JPRB)
  DO JJ=1,1
    DO JF=1,KFC
      PIA(JJ,JF) = (ZNORM/2.0_JPRB) * PIA(JJ,JF)
    ENDDO
  ENDDO
  DO JJ=3,R%NDGL+R%NNOEXTZG+1
    DO JF=1,KFC
      PIA(JJ,JF) = ZNORM * PIA(JJ,JF)
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('ELEINVAD_MOD:ELEINVAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ELEINVAD
END MODULE ELEINVAD_MOD
