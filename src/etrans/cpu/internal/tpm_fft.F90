! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FFT

! Module for Fourier transforms (LAM).
! Provides per-latitude FFT992 tables (TRIGS/NFAX) plus Bluestein fallback
! for lengths incompatible with FFT992, mirroring TPM_FFT992 in trans.

USE PARKIND1, ONLY : JPIM, JPRB
#ifdef WITH_FFT992
USE BLUESTEIN_MOD, ONLY : fftb_type
#endif

IMPLICIT NONE

SAVE

TYPE FFT_TYPE
  REAL(KIND=JPRB),ALLOCATABLE :: TRIGS(:,:)
  INTEGER(KIND=JPIM),ALLOCATABLE :: NFAX(:,:)
  LOGICAL,ALLOCATABLE :: LUSEFFT992(:)
#ifdef WITH_FFT992
  LOGICAL :: LBLUESTEIN=.FALSE.
  TYPE(FFTB_TYPE) :: TB
#endif
END TYPE FFT_TYPE

TYPE(FFT_TYPE),ALLOCATABLE,TARGET :: FFT_RESOL(:)
TYPE(FFT_TYPE),POINTER :: T

END MODULE TPM_FFT
