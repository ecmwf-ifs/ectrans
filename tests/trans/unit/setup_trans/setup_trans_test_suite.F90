! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SETUP_TRANS_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRD
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_END

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_end.h"

! Spectral truncation used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION = 79

! Number of latitudes used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: NDGL = 2 * (TRUNCATION + 1)

CONTAINS

!---------------------------------------------------------------------------------------------------

! Setup fixture
SUBROUTINE SETUP_TEST(LUSE_MPI, NPROC)
  USE UTIL, ONLY: DETECT_MPIRUN

  LOGICAL, INTENT(OUT) :: LUSE_MPI
  INTEGER(KIND=JPIM), INTENT(OUT) :: NPROC

  ! Set up MPI
  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    NPROC = MPL_NPROC()
  ELSE
    NPROC = 1
  ENDIF
END SUBROUTINE SETUP_TEST

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS without first calling SETUP_TRANS0 - should fail
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_WITHOUT_SETUP_TRANS0() RESULT(RET) BIND(C)
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL)
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_WITHOUT_SETUP_TRANS0

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with regular lat-lon grid of 2*(TRUNCATION+1) latitudes
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_BASIC() RESULT(RET) BIND(C)
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)
  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL)
  CALL TRANS_END
  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_BASIC

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with an odd number of latitudes - should fail
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_ODD_NDGL() RESULT(RET) BIND(C)
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)
  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL - 1)
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_ODD_NDGL

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with octahedral grid of 2*(TRUNCATION+1) latitudes 
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_OCTAHEDRAL() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN)
  CALL TRANS_END
  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_OCTAHEDRAL

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with stretch factor passed
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_STRETCHING() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, PSTRET=2.0_JPRD)
  CALL TRANS_END
  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_STRETCHING

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with fast Legendre transform 
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_FLT() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, LDUSEFLT=.TRUE.)
  CALL TRANS_END
  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_FLT

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with all fields passed to FFTW at ocne
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_ALL_FFTW() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, LD_ALL_FFTW=.TRUE.)
  CALL TRANS_END
  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_ALL_FFTW

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with Belusov algorithm
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS_BELUSOV() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I
  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  CALL SETUP_TEST(LUSE_MPI, NPROC)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, LDUSERPNM=.TRUE.)
  CALL TRANS_END
  IF (LUSE_MPI) THEN
    CALL MPL_END
  ENDIF
  RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS_BELUSOV

!---------------------------------------------------------------------------------------------------

END MODULE SETUP_TRANS_TEST_SUITE
