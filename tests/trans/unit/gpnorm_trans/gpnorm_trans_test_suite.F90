! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE GPNORM_TRANS_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB, JPRD
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "dist_grid.h"
#include "gpnorm_trans.h"
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

    ! Set one element to one, all others to zero
    ZGPG(:,:) = 0.0_JPRB
    ZGPG(1,1) = 1.0_JPRB
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

! Test GPNORM_TRANS minimum functionality
INTEGER FUNCTION UNIT_TEST_GPNORM_TRANS_MIN() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:)
  REAL(KIND=JPRB) :: ZAVE(1), ZMIN(1), ZMAX(1)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI

  ! Set up everything
  ! ZGP will contain a single one for one task and zeros elsewhere
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ZGP = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 1)

  ! Make the task a little harder - negate ZGP so there's one negative one and zeros elsewhere
  ZGP = -ZGP

  ! Calculate average, minimum, and maximum of field (LDAVE_ONLY = .FALSE.)
  ! Note: not possible to compute only one of average, min, max - all three must be computed
  ! (unless LDAVE_ONLY = .TRUE. in which case only average is computed)
  CALL GPNORM_TRANS(ZGP, 1, NPROMA, ZAVE, ZMIN, ZMAX, .FALSE.)

  ! If I am task 1, check that the minimum is negative one
  IF (MY_PROC == 1) THEN
    RET = MERGE(0, 1, ZMIN(1) == -1.0_JPRB)
  ELSE
    RET = 0
  ENDIF

  ! Communicate results to other tasks
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_GPNORM_TRANS_MIN

!---------------------------------------------------------------------------------------------------

! Test GPNORM_TRANS maximum functionality
INTEGER FUNCTION UNIT_TEST_GPNORM_TRANS_MAX() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:)
  REAL(KIND=JPRB) :: ZAVE(1), ZMIN(1), ZMAX(1)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI

  ! Set up everything
  ! ZGP will contain a single one for one task and zeros elsewhere
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)
  ZGP = GET_INPUT_FIELD(MY_PROC, NGPTOTG, NGPBLKS, 1)

  ! Calculate average, minimum, and maximum of field (LDAVE_ONLY = .FALSE.)
  ! Note: not possible to compute only one of average, min, max - all three must be computed
  ! (unless LDAVE_ONLY = .TRUE. in which case only average is computed)
  CALL GPNORM_TRANS(ZGP, 1, NPROMA, ZAVE, ZMIN, ZMAX, .FALSE.)

  ! If I am task 1, check that the maximum is one
  IF (MY_PROC == 1) THEN
    RET = MERGE(0, 1, ZMAX(1) == 1.0_JPRB)
  ELSE
    RET = 0
  ENDIF

  ! Communicate results to other tasks
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_GPNORM_TRANS_MAX

!---------------------------------------------------------------------------------------------------

! Test GPNORM_TRANS average functionality
INTEGER FUNCTION UNIT_TEST_GPNORM_TRANS_AVE() RESULT(RET) BIND(C)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:)
  REAL(KIND=JPRB) :: ZAVE(1), ZMIN(1), ZMAX(1)
  INTEGER(KIND=JPIM) :: NSPEC2, NGPTOTG, NGPTOT, MY_PROC, NGPBLKS
  LOGICAL :: LUSE_MPI

  ! Set up everything
  CALL SETUP_TEST(NSPEC2, NGPTOTG, NGPTOT, NGPBLKS, LUSE_MPI, MY_PROC)

  ! Initialise distributed fields
  ALLOCATE(ZGP(NPROMA,1,NGPBLKS))

  ! Set all values to one
  ZGP(:,:,:) = 1.0_JPRB

  ! Calculate average, minimum, and maximum of field (LDAVE_ONLY = .FALSE.)
  ! Note: not possible to compute only one of average, min, max - all three must be computed
  ! (unless LDAVE_ONLY = .TRUE. in which case only average is computed)
  CALL GPNORM_TRANS(ZGP, 1, NPROMA, ZAVE, ZMIN, ZMAX, .FALSE.)

  ! If I am task 1, check that the average is one
  IF (MY_PROC == 1) THEN
    RET = MERGE(0, 1, APPROX_EQ(ZAVE(1), 1.0_JPRB))
  ELSE
    RET = 0
  ENDIF

  ! Communicate results to other tasks
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(ZGP)
  CALL CLEANUP_TEST(LUSE_MPI)
END FUNCTION UNIT_TEST_GPNORM_TRANS_AVE

!---------------------------------------------------------------------------------------------------

END MODULE GPNORM_TRANS_TEST_SUITE
