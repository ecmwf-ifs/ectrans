! (C) Copyright 2024- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! ecTrans demonstrator program
! see https://sites.ecmwf.int/docs/ectrans/page/usage.html for a full explanation

PROGRAM ECTRANS_DEMONSTRATOR

! Import 4-byte float and 4-byte integer types
USE PARKIND1, ONLY: JPRM, JPIM

IMPLICIT NONE

! Include header files for every ecTrans subroutine we want to call
#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"

! Spectral truncation
INTEGER, PARAMETER :: TRUNC = 79

! Arrays for storing our field in spectral space and grid point space
REAL(KIND=JPRM), ALLOCATABLE :: SPECTRAL_FIELD(:,:)
REAL(KIND=JPRM), ALLOCATABLE :: SPECTRAL_FIELD_2(:,:)
REAL(KIND=JPRM), ALLOCATABLE :: GRID_POINT_FIELD(:,:,:)

! Spectral indices array for initialising spectral field
INTEGER(KIND=JPIM) :: SPECTRAL_INDICES(0:TRUNC)

! Dimensions of our arrays in spectral space and grid point space
INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NGPTOT

INTEGER(KIND=JPIM) :: M, N, I

! Initialise ecTrans (resolution-agnostic aspects)
CALL SETUP_TRANS0(LDMPOFF=.TRUE.)

! Initialise ecTrans (resolution-specific aspects)
CALL SETUP_TRANS(KSMAX=TRUNC, KDGL=2 * (TRUNC + 1))

! Inquire about the dimensions in spectral space and grid point space
CALL TRANS_INQ(KSPEC2=NSPEC2, KGPTOT=NGPTOT)

WRITE(*,*) "NSPEC2 = ", NSPEC2
WRITE(*,*) "NGPTOT = ", NGPTOT

! Allocate our work arrays
ALLOCATE(SPECTRAL_FIELD(1,NSPEC2))
ALLOCATE(SPECTRAL_FIELD_2(1,NSPEC2))
ALLOCATE(GRID_POINT_FIELD(NGPTOT,1,1))

! Initialise our spectral field array
CALL TRANS_INQ(KASM0=SPECTRAL_INDICES)
M = 3; N = 8
I = SPECTRAL_INDICES(M) + 2 * (N - M) + 1
SPECTRAL_FIELD(:,:) = 0.0
SPECTRAL_FIELD(1,I) = 1.0

! Perform an inverse transform
CALL INV_TRANS(PSPSCALAR=SPECTRAL_FIELD, PGP=GRID_POINT_FIELD)

! Dump the grid point field to an unformatted binary file
OPEN(7, FILE="grid_point_field.dat", FORM="UNFORMATTED")
WRITE(7) GRID_POINT_FIELD(:,1,1)
CLOSE(7)

! Perform a direct transform back
CALL DIR_TRANS(PGP=GRID_POINT_FIELD, PSPSCALAR=SPECTRAL_FIELD_2)

! Compute error between first and second spectral field arrays
WRITE(*,*) "Error = ", NORM2(SPECTRAL_FIELD_2 - SPECTRAL_FIELD)

END PROGRAM ECTRANS_DEMONSTRATOR
