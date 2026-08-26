! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM TEST_MULTIPLE_PRECISIONS

USE PARKIND1, ONLY: JPRM, JPRD, JPIM

IMPLICIT NONE

! This one is "precision agnostic" -> there is no "single-precision version of SETUP_TRANS0"
#include "setup_trans0.h"

! The following ARE dependent on precision
! By including _sp or _dp versions, we can have multiple precisions in the same program
#include "setup_trans_sp.h"
#include "setup_trans_dp.h"
#include "trans_inq_sp.h"
#include "inv_trans_sp.h"
#include "inv_trans_dp.h"
#include "trans_end_sp.h"
#include "trans_end_dp.h"

! STDOUT
INTEGER, PARAMETER :: NOUT = 6

! Spectral truncation
INTEGER, PARAMETER :: TRUNC = 79

REAL(KIND=JPRM), ALLOCATABLE :: SPECTRAL_FIELD_SP(:,:)
REAL(KIND=JPRM), ALLOCATABLE :: GRID_POINT_FIELD_SP(:,:,:)
REAL(KIND=JPRD), ALLOCATABLE :: SPECTRAL_FIELD_DP(:,:)
REAL(KIND=JPRD), ALLOCATABLE :: GRID_POINT_FIELD_DP(:,:,:)

! Dimensions of our arrays in spectral space and grid point space
INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NGPTOT

INTEGER(KIND=JPIM) :: IRESOL_SP, IRESOL_DP

! This call only needs to happen once
CALL SETUP_TRANS0(KOUT=NOUT, LDMPOFF=.TRUE., KMAX_RESOL=2)

! Each precision needs a separate call to SETUP_TRANS
CALL SETUP_TRANS_SP(KSMAX=TRUNC, KDGL=2 * (TRUNC + 1), KRESOL=IRESOL_SP)
WRITE(NOUT,*) "Single-precision handle = ", IRESOL_SP
CALL SETUP_TRANS_DP(KSMAX=TRUNC, KDGL=2 * (TRUNC + 1), KRESOL=IRESOL_DP)
WRITE(NOUT,*) "Double-precision handle = ", IRESOL_DP

! Inquire about the dimensions in spectral space and grid point space
CALL TRANS_INQ_SP(KSPEC2=NSPEC2, KGPTOT=NGPTOT)

! Allocate our work arrays
ALLOCATE(SPECTRAL_FIELD_SP(1,NSPEC2))
ALLOCATE(GRID_POINT_FIELD_SP(NGPTOT,1,1))
ALLOCATE(SPECTRAL_FIELD_DP(1,NSPEC2))
ALLOCATE(GRID_POINT_FIELD_DP(NGPTOT,1,1))

! Initialise our spectral field arrays
SPECTRAL_FIELD_SP(:,:) = 1.0_JPRM
SPECTRAL_FIELD_DP(:,:) = 1.0_JPRD

! Perform an inverse transform
CALL INV_TRANS_SP(PSPSCALAR=SPECTRAL_FIELD_SP, PGP=GRID_POINT_FIELD_SP)
CALL INV_TRANS_DP(PSPSCALAR=SPECTRAL_FIELD_DP, PGP=GRID_POINT_FIELD_DP)

WRITE(NOUT,*) "Grid point norm (single precision) = ", SQRT(SUM(GRID_POINT_FIELD_SP(:,:,:)**2.0_JPRM))
WRITE(NOUT,*) "Grid point norm (double precision) = ", SQRT(SUM(GRID_POINT_FIELD_DP(:,:,:)**2.0_JPRD))

CALL TRANS_END_SP
CALL TRANS_END_DP

WRITE(NOUT,*) "Finished"

END PROGRAM TEST_MULTIPLE_PRECISIONS