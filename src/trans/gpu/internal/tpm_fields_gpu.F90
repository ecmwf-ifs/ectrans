! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2024- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FIELDS_GPU

USE EC_PARKIND, ONLY: JPRD, JPRBT

IMPLICIT NONE

SAVE

TYPE FIELDS_GPU_TYPE
! scratch arrays for ltinv and ltdir and associated dimension variables
REAL(KIND=JPRBT),ALLOCATABLE :: ZAA(:,:,:)  !! JPRL for 1/2
REAL(KIND=JPRBT),ALLOCATABLE :: ZAS(:,:,:)  !! JPRL for 1/2

! for m=0 in ledir_mod:
REAL(KIND=JPRD),ALLOCATABLE :: ZAA0(:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAS0(:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZEPSNM(:,:)
END TYPE FIELDS_GPU_TYPE

TYPE(FIELDS_GPU_TYPE),ALLOCATABLE,TARGET :: FIELDS_GPU_RESOL(:)
TYPE(FIELDS_GPU_TYPE),POINTER     :: FG

END MODULE TPM_FIELDS_GPU
