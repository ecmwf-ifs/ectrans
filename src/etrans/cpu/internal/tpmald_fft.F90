! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE TPMALD_FFT

! Module for Fourier transforms.

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE ALDFFT_TYPE
REAL(KIND=JPRB)   ,POINTER :: TRIGSE(:) ! list of trigonometric function values
INTEGER(KIND=JPIM),POINTER :: NFAXE(:)  ! list of factors of truncation
LOGICAL                    :: LFFT992=.FALSE.
END TYPE ALDFFT_TYPE

TYPE(ALDFFT_TYPE),ALLOCATABLE,TARGET :: ALDFFT_RESOL(:)
TYPE(ALDFFT_TYPE),POINTER     :: TALD

END MODULE TPMALD_FFT
