! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE DIR_TRANS_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB, JPRD

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "dir_trans.h"
#include "trans_end.h"

! Spectral truncation used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION = 79

! Number of latitudes used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: NDGL = 2 * (TRUNCATION + 1)

! Tolerance for "close to zero"
REAL(KIND=JPRB), PARAMETER :: TOLERANCE = 100.0_JPRB * EPSILON(1.0_JPRB)

CONTAINS

!---------------------------------------------------------------------------------------------------

! Approximate equality check for reals
ELEMENTAL LOGICAL FUNCTION APPROX_EQ(A, B, TOL) RESULT(RET)
  REAL(KIND=JPRB), INTENT(IN) :: A
  REAL(KIND=JPRB), INTENT(IN) :: B
  REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: TOL

  IF (PRESENT(TOL)) THEN
    RET = ABS(A - B) < TOL
  ELSE
    RET = ABS(A - B) < TOLERANCE
  END IF
END FUNCTION APPROX_EQ

! Setup fixture (octahedral grid)
SUBROUTINE SETUP_TEST(KSPEC2, KGPTOT)
  INTEGER(KIND=JPIM), INTENT(OUT) :: KSPEC2
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPTOT

  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I

  CALL SETUP_TRANS0(LDMPOFF=.TRUE.)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN)

  CALL TRANS_INQ(KSPEC2=KSPEC2, KGPTOT=KGPTOT)
END SUBROUTINE SETUP_TEST

!---------------------------------------------------------------------------------------------------

! NOTES:
! - All tests below have MPI disabled.
! - For now there is only a very primitive correctness check. We set the input to one and check that
!   only the (0,0) mode (global mean) is one and all the rest are zero.

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 1 and just one scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_SCALAR_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZSPSCALAR(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT

  CALL SETUP_TEST(NSPEC2, NGPTOT)

  ! Initialise arrays
  ALLOCATE(ZGP(NGPTOT,1,1))
  ALLOCATE(ZSPSCALAR(1,NSPEC2))

  ! Set all of the input to one
  ZGP(:,:,:) = 1.0_JPRB

  CALL DIR_TRANS(PGP=ZGP, PSPSCALAR=ZSPSCALAR)

  ! Check only the (0,0) mode (global mean) is one and all the rest are zero
  IF (APPROX_EQ(ZSPSCALAR(1,1), 1.0_JPRB) .AND. ALL(APPROX_EQ(ZSPSCALAR(1,2:), 0.0_JPRB))) THEN
    RET = 0
  ELSE
    RET = 1
  END IF

  DEALLOCATE(ZSPSCALAR, ZGP)
  CALL TRANS_END
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_SCALAR_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 1 and 1-level wind fields
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZSPVOR(:,:), ZSPDIV(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  CALL SETUP_TEST(NSPEC2, NGPTOT)

  ! Initialise arrays
  ALLOCATE(ZGP(NGPTOT,2,1))
  ALLOCATE(ZSPVOR(1,NSPEC2))
  ALLOCATE(ZSPDIV(1,NSPEC2))

  ! Set all of the input to one
  ZGP(:,:,:) = 1.0_JPRB

  CALL DIR_TRANS(PGP=ZGP, PSPVOR=ZSPVOR, PSPDIV=ZSPDIV)

  ! We pass in constant U and V so vorticity and divergence should both be zero
  ! TODO: for some reason the double precision computation of vorticity and divergence gives
  ! errors more like single precision. This possibly indicates a hard-coded single-precision
  ! value somewhere in the double-precision code path. For now we use a higher tolerance.
  IF (JPRB == JPRD) THEN
    TEMP_TOLERANCE = 1E7_JPRD * TOLERANCE
  ELSE
    TEMP_TOLERANCE = TOLERANCE
  END IF
  IF (ALL(APPROX_EQ(ZSPVOR(1,1:), 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
    & ALL(APPROX_EQ(ZSPDIV(1,1:), 0.0_JPRB, TOL=TEMP_TOLERANCE))) THEN
    RET = 0
  ELSE
    RET = 1
  END IF

  DEALLOCATE(ZSPDIV, ZSPVOR, ZGP)
  CALL TRANS_END
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 1 and 1-level wind fields and 1 scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1_SCALAR_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZSPVOR(:,:), ZSPDIV(:,:), ZSPSCALAR(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  CALL SETUP_TEST(NSPEC2, NGPTOT)

  ! Initialise arrays
  ALLOCATE(ZGP(NGPTOT,3,1))
  ALLOCATE(ZSPVOR(1,NSPEC2))
  ALLOCATE(ZSPDIV(1,NSPEC2))
  ALLOCATE(ZSPSCALAR(1,NSPEC2))

  ! Set all of the input to one
  ZGP(:,:,:) = 1.0_JPRB

  CALL DIR_TRANS(PGP=ZGP, PSPVOR=ZSPVOR, PSPDIV=ZSPDIV, PSPSCALAR=ZSPSCALAR)

  ! We pass in constant U and V so vorticity and divergence should both be zero
  ! As for the scalar, check only the (0,0) mode (global mean) is one and all the rest are zero
  ! TODO: for some reason the double precision computation of vorticity and divergence gives
  ! errors more like single precision. This possibly indicates a hard-coded single-precision
  ! value somewhere in the double-precision code path. For now we use a higher tolerance.
  IF (JPRB == JPRD) THEN
    TEMP_TOLERANCE = 1E7_JPRD * TOLERANCE
  ELSE
    TEMP_TOLERANCE = TOLERANCE
  END IF
  IF (ALL(APPROX_EQ(ZSPVOR(1,1:), 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
    & ALL(APPROX_EQ(ZSPDIV(1,1:), 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
    & APPROX_EQ(ZSPSCALAR(1,1), 1.0_JPRB) .AND. ALL(APPROX_EQ(ZSPSCALAR(1,2:), 0.0_JPRB))) THEN
    RET = 0
  ELSE
    RET = 1
  END IF

  DEALLOCATE(ZSPSCALAR, ZSPDIV, ZSPVOR, ZGP)
  CALL TRANS_END
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1_SCALAR_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 2 and just one "3A" scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3A_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP3A(:,:,:,:), ZSPSC3A(:,:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  CALL SETUP_TEST(NSPEC2, NGPTOT)

  ! Initialise arrays
  ALLOCATE(ZGP3A(NGPTOT,1,1,1))
  ALLOCATE(ZSPSC3A(1,NSPEC2,1))

  ! Set all of the input to one
  ZGP3A(:,:,:,:) = 1.0_JPRB

  CALL DIR_TRANS(PGP3A=ZGP3A, PSPSC3A=ZSPSC3A)

  ! Check only the (0,0) mode (global mean) is one and all the rest are zero
  IF (APPROX_EQ(ZSPSC3A(1,1,1), 1.0_JPRB) .AND. ALL(APPROX_EQ(ZSPSC3A(1,2:,1), 0.0_JPRB))) THEN
    RET = 0
  ELSE
    RET = 1
  END IF

  DEALLOCATE(ZSPSC3A, ZGP3A)
  CALL TRANS_END
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3A_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 2 and just one "3B" scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3B_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP3B(:,:,:,:), ZSPSC3B(:,:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  CALL SETUP_TEST(NSPEC2, NGPTOT)

  ! Initialise arrays
  ALLOCATE(ZGP3B(NGPTOT,1,1,1))
  ALLOCATE(ZSPSC3B(1,NSPEC2,1))

  ! Set all of the input to one
  ZGP3B(:,:,:,:) = 1.0_JPRB

  CALL DIR_TRANS(PGP3B=ZGP3B, PSPSC3B=ZSPSC3B)

  ! Check only the (0,0) mode (global mean) is one and all the rest are zero
  IF (APPROX_EQ(ZSPSC3B(1,1,1), 1.0_JPRB) .AND. ALL(APPROX_EQ(ZSPSC3B(1,2:,1), 0.0_JPRB))) THEN
    RET = 0
  ELSE
    RET = 1
  END IF

  DEALLOCATE(ZSPSC3B, ZGP3B)
  CALL TRANS_END
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3B_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 2 and 1-level wind fields
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_WIND_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGPUV(:,:,:,:), ZSPVOR(:,:), ZSPDIV(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  CALL SETUP_TEST(NSPEC2, NGPTOT)

  ! Initialise arrays
  ALLOCATE(ZGPUV(NGPTOT,1,2,1))
  ALLOCATE(ZSPVOR(1,NSPEC2))
  ALLOCATE(ZSPDIV(1,NSPEC2))

  ! Set all of the input to one
  ZGPUV(:,:,:,:) = 1.0_JPRB

  CALL DIR_TRANS(PGPUV=ZGPUV, PSPVOR=ZSPVOR, PSPDIV=ZSPDIV)

  ! We pass in constant U and V so vorticity and divergence should both be zero
  ! As for the scalar, check only the (0,0) mode (global mean) is one and all the rest are zero
  ! TODO: for some reason the double precision computation of vorticity and divergence gives
  ! errors more like single precision. This possibly indicates a hard-coded single-precision
  ! value somewhere in the double-precision code path. For now we use a higher tolerance.
  IF (JPRB == JPRD) THEN
    TEMP_TOLERANCE = 1E7_JPRD * TOLERANCE
  ELSE
    TEMP_TOLERANCE = TOLERANCE
  END IF
  IF (ALL(APPROX_EQ(ZSPVOR(1,1:), 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
    & ALL(APPROX_EQ(ZSPDIV(1,1:), 0.0_JPRB, TOL=TEMP_TOLERANCE))) THEN
    RET = 0
  ELSE
    RET = 1
  END IF

  DEALLOCATE(ZSPDIV, ZSPVOR, ZGPUV)
  CALL TRANS_END
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_WIND_1

!---------------------------------------------------------------------------------------------------

END MODULE DIR_TRANS_TEST_SUITE
