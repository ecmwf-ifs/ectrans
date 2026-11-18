! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FFT992

USE PARKIND1, ONLY : JPIM, JPRB
#ifdef WITH_FFT992
USE BLUESTEIN_MOD, ONLY : fftb_type
#endif

IMPLICIT NONE

SAVE

TYPE FFT992_TYPE
  REAL(KIND=JPRB),ALLOCATABLE :: TRIGS(:,:)
  INTEGER(KIND=JPIM),ALLOCATABLE :: NFAX(:,:)
  LOGICAL,ALLOCATABLE :: LUSEFFT992(:)
#ifdef WITH_FFT992
  LOGICAL :: LBLUESTEIN=.FALSE.
#endif
  LOGICAL :: LFFT992=.FALSE.
#ifdef WITH_FFT992
  TYPE(FFTB_TYPE) :: TB
#endif
END TYPE FFT992_TYPE

TYPE(FFT992_TYPE),ALLOCATABLE,TARGET :: FFT992_RESOL(:)
TYPE(FFT992_TYPE),POINTER :: T992

END MODULE TPM_FFT992
