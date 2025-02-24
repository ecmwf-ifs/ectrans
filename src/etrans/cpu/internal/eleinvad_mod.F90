! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


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
USE TPM_FFTW     ,ONLY : TW, EXEC_EFFTW
USE TPMALD_DIM      ,ONLY : RALD
#ifdef WITH_FFT992
USE TPMALD_FFT      ,ONLY : TALD
#endif
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
#ifdef WITH_FFT992
  IF( TALD%LFFT992 )THEN
    CALL FFT992(PIA,TALD%TRIGSE,TALD%NFAXE,1,ICLEN,IRLEN,KFC,ITYPE)
  ELSE
#endif
    IOFF=1
    CALL EXEC_EFFTW(ITYPE,IRLEN,ICLEN,IOFF,KFC,TW%LALL_FFTW,PIA)
#ifdef WITH_FFT992
  ENDIF
#endif
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
