! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ELEDIR_MOD
CONTAINS
SUBROUTINE ELEDIR(KM,KFC,KLED2,PFFT)

!**** *ELEDIR* - Direct meridional transform.

!     Purpose.
!     --------
!        Direct meridional tranform of state variables.

!**   Interface.
!     ----------
!        CALL ELEDIR(...)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFC - number of field to transform
!                              PAIA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSIA - symmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              POA1 -  spectral
!                              fields for zonal wavenumber KM
!                              PLEPO - Legendre polonomials

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - matrix multiply
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - NTMAX instead of NSMAX
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        R. El Khatib 01-Sep-2015 support for FFTW transforms
!        R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_DIM         ,ONLY : R
!USE TPM_GEOMETRY
!USE TPM_TRANS
#ifdef WITH_FFT992
USE TPMALD_FFT      ,ONLY : TALD
#endif
USE TPM_FFTW     ,ONLY : TW, EXEC_EFFTW
USE TPMALD_DIM      ,ONLY : RALD
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM,KFC,KLED2
REAL(KIND=JPRB) ,   INTENT(INOUT)  :: PFFT(:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN, IOFF, ITYPE
!     ------------------------------------------------------------------

!*       1.       PERFORM FOURIER TRANFORM.
!                 --------------------------

IF (KFC>0) THEN
  ITYPE=-1
  IRLEN=R%NDGL+R%NNOEXTZG
  ICLEN=RALD%NDGLSUR+R%NNOEXTZG
#ifdef WITH_FFT992
  IF( TALD%LFFT992 )THEN
    CALL FFT992(PFFT,TALD%TRIGSE,TALD%NFAXE,1,ICLEN,IRLEN,KFC,ITYPE)
#endif
    IOFF=1
    CALL EXEC_EFFTW(ITYPE,IRLEN,ICLEN,IOFF,KFC,TW%LALL_FFTW,PFFT)
#ifdef WITH_FFT992
  ENDIF
#endif
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE ELEDIR
END MODULE ELEDIR_MOD
