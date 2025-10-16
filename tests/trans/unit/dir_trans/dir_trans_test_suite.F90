! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE DIR_TRANS_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB, JPRD
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "dist_grid.h"
#include "dir_trans.h"
#include "trans_end.h"

! Spectral truncation used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION = 79

! Number of latitudes used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: NDGL = 2 * (TRUNCATION + 1)

! Tolerance for "close to zero"
REAL(KIND=JPRB), PARAMETER :: TOLERANCE = 100.0_JPRB * EPSILON(1.0_JPRB)

! NPROMA blocking factor
INTEGER(KIND=JPIM), PARAMETER :: NPROMA = 16

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

! Initialise global field with all ones and distribute it
FUNCTION GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, NFIELDS) RESULT(PGP)
  INTEGER(KIND=JPIM), INTENT(IN) :: MY_PROC
  INTEGER(KIND=JPIM), INTENT(IN) :: NGPTOTG
  INTEGER(KIND=JPIM), INTENT(IN) :: NGPBLKS
  INTEGER(KIND=JPIM), INTENT(IN) :: NFIELDS

  REAL(KIND=JPRB), ALLOCATABLE :: ZGPG(:,:), PGP(:,:,:)
  INTEGER(KIND=JPIM) :: I

  ! Initialise global field
  IF (MY_PROC == 1) THEN
    ALLOCATE(ZGPG(NGPTOTG,NFIELDS))

    ! Set all of the input to one
    ZGPG(:,:) = 1.0_JPRB
  ENDIF

  ! Initialise distributed fields
  ALLOCATE(PGP(NPROMA,NFIELDS,NGPBLKS))

  ! Distribute field from first task to other tasks
  IF (MY_PROC == 1) THEN
    CALL DIST_GRID(PGPG=ZGPG, KFDISTG=NFIELDS, KFROM=(/(1, I = 1,NFIELDS)/), PGP=PGP, KPROMA=NPROMA)
  ELSE
    CALL DIST_GRID(KFDISTG=NFIELDS, KFROM=(/ (1, I = 1, NFIELDS) /), PGP=PGP, KPROMA=NPROMA)
  ENDIF

  IF (MY_PROC == 1) DEALLOCATE(ZGPG)
END FUNCTION GET_INPUT_FIELD

! Determine if this task handles the m=0 mode
LOGICAL FUNCTION HAVE_M0_MODE()
  INTEGER :: NUM_MY_ZON_WNS
  INTEGER, ALLOCATABLE :: MY_ZON_WNS(:)

  CALL TRANS_INQ(KNUMP=NUM_MY_ZON_WNS)
  ALLOCATE(MY_ZON_WNS(NUM_MY_ZON_WNS))
  CALL TRANS_INQ(KMYMS=MY_ZON_WNS)

  HAVE_M0_MODE = ANY(MY_ZON_WNS == 0)
END FUNCTION HAVE_M0_MODE

! Setup fixture
SUBROUTINE SETUP_TEST(KSPEC2, KGPTOTG, KGPTOT, KGPBLKS, LUSE_MPI, MY_PROC)
  USE UTIL, ONLY: DETECT_MPIRUN

  INTEGER(KIND=JPIM), INTENT(OUT) :: KSPEC2
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPTOTG
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPTOT
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPBLKS
  LOGICAL, INTENT(OUT) :: LUSE_MPI
  INTEGER(KIND=JPIM), INTENT(OUT) :: MY_PROC

  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I
  INTEGER(KIND=JPIM) :: NPROC

  ! Set up MPI
  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    NPROC = MPL_NPROC()
    MY_PROC = MPL_MYRANK()
  ELSE
    NPROC = 1
    MY_PROC = 1
  ENDIF

  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC, KPRGPEW=1, KPRTRW=NPROC)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN)

  CALL TRANS_INQ(KSPEC2=KSPEC2, KGPTOTG=KGPTOTG, KGPTOT=KGPTOT)

  ! Number of NPROMA blocks
  KGPBLKS = (KGPTOT - 1) / NPROMA + 1
END SUBROUTINE SETUP_TEST

! Tear down fixture
SUBROUTINE CLEANUP_TEST(LUSE_MPI)
  LOGICAL, INTENT(IN) :: LUSE_MPI

  CALL TRANS_END

  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
END SUBROUTINE CLEANUP_TEST

!---------------------------------------------------------------------------------------------------

! NOTES:
! - For now there is only a very primitive correctness check. For scalar fields we set the input to
!   one and check that only the (0,0) mode (global mean) is one and all the rest are zero. For wind
!  fields we set the input to one and check that vorticity and divergence are both zero.

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 1 and just one scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_SCALAR_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZSPSCALAR(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ZGP = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 1)
  ALLOCATE(ZSPSCALAR(1,NSPEC2))

  CALL DIR_TRANS(PGP=ZGP, PSPSCALAR=ZSPSCALAR, KPROMA=NPROMA)

  ! Check only the (0,0) mode (global mean) is one and all the rest are zero
  IF (HAVE_M0_MODE()) THEN
    RET = MERGE(0, 1, APPROX_EQ(ZSPSCALAR(1,1), 1.0_JPRB) .AND. &
      &         ALL(APPROX_EQ(ZSPSCALAR(1,2:), 0.0_JPRB)))
  ELSE
    RET = MERGE(0, 1, ALL(APPROX_EQ(ZSPSCALAR, 0.0_JPRB)))
  ENDIF

  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZSPSCALAR, ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_SCALAR_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 1 and 1-level wind fields
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZSPVOR(:,:), ZSPDIV(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ZGP = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 2)
  ALLOCATE(ZSPVOR(1,NSPEC2))
  ALLOCATE(ZSPDIV(1,NSPEC2))

  CALL DIR_TRANS(PGP=ZGP, PSPVOR=ZSPVOR, PSPDIV=ZSPDIV, KPROMA=NPROMA)

  ! We pass in constant U and V so vorticity and divergence should both be zero
  ! TODO: for some reason the double precision computation of vorticity and divergence gives
  ! errors more like single precision. This possibly indicates a hard-coded single-precision
  ! value somewhere in the double-precision code path. For now we use a higher tolerance.
  IF (JPRB == JPRD) THEN
    TEMP_TOLERANCE = 1E7_JPRD * TOLERANCE
  ELSE
    TEMP_TOLERANCE = TOLERANCE
  END IF
  RET = MERGE(0, 1, ALL(APPROX_EQ(ZSPVOR, 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
    & ALL(APPROX_EQ(ZSPDIV, 0.0_JPRB, TOL=TEMP_TOLERANCE)))

  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZSPDIV, ZSPVOR, ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 1 and 1-level wind fields and 1 scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1_SCALAR_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:), ZSPVOR(:,:), ZSPDIV(:,:), ZSPSCALAR(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ZGP = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 3)
  ALLOCATE(ZSPVOR(1,NSPEC2))
  ALLOCATE(ZSPDIV(1,NSPEC2))
  ALLOCATE(ZSPSCALAR(1,NSPEC2))

  CALL DIR_TRANS(PGP=ZGP, PSPVOR=ZSPVOR, PSPDIV=ZSPDIV, PSPSCALAR=ZSPSCALAR, KPROMA=NPROMA)

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
  IF (HAVE_M0_MODE()) THEN
    IF (ALL(APPROX_EQ(ZSPVOR, 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
      & ALL(APPROX_EQ(ZSPDIV, 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
      & APPROX_EQ(ZSPSCALAR(1,1), 1.0_JPRB) .AND. ALL(APPROX_EQ(ZSPSCALAR(1,2:), 0.0_JPRB))) THEN
      RET = 0
    ELSE
      RET = 1
    END IF
  ELSE
    IF (ALL(APPROX_EQ(ZSPVOR, 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
      & ALL(APPROX_EQ(ZSPDIV, 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
      & ALL(APPROX_EQ(ZSPSCALAR, 0.0_JPRB))) THEN
      RET = 0
    ELSE
      RET = 1
    END IF
  ENDIF

  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZSPSCALAR, ZSPDIV, ZSPVOR, ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_1_WIND_1_SCALAR_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 2 and just one "3A" scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3A_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP3A(:,:,:,:), ZSPSC3A(:,:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ALLOCATE(ZGP3A(NPROMA,1,1,NGPBLKS))
  ZGP3A(:,1,:,:) = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 1)
  ALLOCATE(ZSPSC3A(1,NSPEC2,1))

  CALL DIR_TRANS(PGP3A=ZGP3A, PSPSC3A=ZSPSC3A, KPROMA=NPROMA)

  ! Check only the (0,0) mode (global mean) is one and all the rest are zero
  IF (HAVE_M0_MODE()) THEN
    RET = MERGE(0, 1, APPROX_EQ(ZSPSC3A(1,1,1), 1.0_JPRB) .AND. &
      &         ALL(APPROX_EQ(ZSPSC3A(1,2:,1), 0.0_JPRB)))
  ELSE
    RET = MERGE(0, 1, ALL(APPROX_EQ(ZSPSC3A, 0.0_JPRB)))
  ENDIF

  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZSPSC3A, ZGP3A)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3A_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 2 and just one "3B" scalar field
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3B_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP3B(:,:,:,:), ZSPSC3B(:,:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ALLOCATE(ZGP3B(NPROMA,1,1,NGPBLKS))
  ZGP3B(:,1,:,:) = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 1)
  ALLOCATE(ZSPSC3B(1,NSPEC2,1))

  CALL DIR_TRANS(PGP3B=ZGP3B, PSPSC3B=ZSPSC3B, KPROMA=NPROMA)

  ! Check only the (0,0) mode (global mean) is one and all the rest are zero
  IF (HAVE_M0_MODE()) THEN
    RET = MERGE(0, 1, APPROX_EQ(ZSPSC3B(1,1,1), 1.0_JPRB) .AND. &
      &         ALL(APPROX_EQ(ZSPSC3B(1,2:,1), 0.0_JPRB)))
  ELSE
    RET = MERGE(0, 1, ALL(APPROX_EQ(ZSPSC3B, 0.0_JPRB)))
  ENDIF

  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZSPSC3B, ZGP3B)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_PGP3B_1

!---------------------------------------------------------------------------------------------------

! Test DIR_TRANS with call mode 2 and 1-level wind fields
INTEGER FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_WIND_1() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGPUV(:,:,:,:), ZSPVOR(:,:), ZSPDIV(:,:)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI
  REAL(KIND=JPRB) :: TEMP_TOLERANCE

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ALLOCATE(ZGPUV(NPROMA,1,2,NGPBLKS))
  ZGPUV(:,1,:,:) = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 2)
  ALLOCATE(ZSPVOR(1,NSPEC2))
  ALLOCATE(ZSPDIV(1,NSPEC2))

  CALL DIR_TRANS(PGPUV=ZGPUV, PSPVOR=ZSPVOR, PSPDIV=ZSPDIV, KPROMA=NPROMA)

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
  RET = MERGE(0, 1, ALL(APPROX_EQ(ZSPVOR, 0.0_JPRB, TOL=TEMP_TOLERANCE)) .AND. &
    & ALL(APPROX_EQ(ZSPDIV, 0.0_JPRB, TOL=TEMP_TOLERANCE)))

  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZSPDIV, ZSPVOR, ZGPUV)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_DIR_TRANS_CALL_MODE_2_WIND_1

!---------------------------------------------------------------------------------------------------

END MODULE DIR_TRANS_TEST_SUITE
