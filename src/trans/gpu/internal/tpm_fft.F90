! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FFT
USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT

! Module for Fourier transforms.

IMPLICIT NONE

SAVE

TYPE FFT_TYPE
  REAL(KIND=JPRBT)  ,ALLOCATABLE :: TRIGS(:,:) ! list of trigonometric function values
  INTEGER(KIND=JPIM),ALLOCATABLE :: NFAX(:,:)  ! list of factors of truncation
END TYPE FFT_TYPE

TYPE(FFT_TYPE),ALLOCATABLE,TARGET :: FFT_RESOL(:)
TYPE(FFT_TYPE),POINTER     :: T


END MODULE TPM_FFT
