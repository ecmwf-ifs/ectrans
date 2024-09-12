! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FIELDS_FLAT

USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT, JPRD

IMPLICIT NONE

SAVE

! flat copies of the fields defined in TPM_FIELDS
REAL(KIND=JPRD)  ,ALLOCATABLE :: F_RW(:)     ! Weights of the Gaussian quadrature
REAL(KIND=JPRBT) ,ALLOCATABLE :: F_RLAPIN(:) ! eigen-values of the inverse Laplace operator
REAL(KIND=JPRD)  ,ALLOCATABLE :: F_RACTHE(:) ! eigen-values of the inverse Laplace operator

! scratch arrays for ltinv and ltdir and associated dimension variables

REAL(KIND=JPRBT),ALLOCATABLE :: ZAA(:,:,:)  !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE :: ZAS(:,:,:)  !! JPRL for 1/2

! for m=0 in ledir_mod:
REAL(KIND=JPRD),ALLOCATABLE :: ZAA0(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAS0(:,:)
INTEGER(KIND=JPIM) :: KMLOC0

REAL(KIND=JPRBT),ALLOCATABLE :: ZEPSNM(:,:)

END MODULE TPM_FIELDS_FLAT
