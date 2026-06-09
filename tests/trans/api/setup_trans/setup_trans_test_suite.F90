! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SETUP_TRANS_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRD

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_end.h"

! Spectral truncation used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION = 79

! Number of latitudes used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: NDGL = 2 * (TRUNCATION + 1)

LOGICAL :: LUSE_MPI
INTEGER(KIND=JPIM) :: NPROC

CONTAINS

!---------------------------------------------------------------------------------------------------

! Setup fixture
SUBROUTINE SETUP_TEST(LCALL_SETUP_TRANS0)
  USE UTIL, ONLY: DETECT_MPIRUN
  USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC

  LOGICAL, INTENT(IN), OPTIONAL :: LCALL_SETUP_TRANS0

  LOGICAL :: LLCALL_SETUP_TRANS0 = .TRUE.

  ! Set up MPI
  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    NPROC = MPL_NPROC()
  ELSE
    NPROC = 1
  ENDIF

  IF (PRESENT(LCALL_SETUP_TRANS0)) THEN
    LLCALL_SETUP_TRANS0 = LCALL_SETUP_TRANS0
  END IF
  IF (LLCALL_SETUP_TRANS0) CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC)
END SUBROUTINE SETUP_TEST

!---------------------------------------------------------------------------------------------------

! Teardown fixture
SUBROUTINE END_TEST
  USE MPL_MODULE, ONLY: MPL_END

  CALL TRANS_END
  IF (LUSE_MPI) THEN
    CALL MPL_END(LDMEMINFO=.FALSE.)
  ENDIF
END SUBROUTINE END_TEST

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS without first calling SETUP_TRANS0 - should fail
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_WITHOUT_SETUP_TRANS0() RESULT(RET) BIND(C)

  CALL SETUP_TEST(LCALL_SETUP_TRANS0=.FALSE.)

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL)

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_WITHOUT_SETUP_TRANS0

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with regular lat-lon grid of 2*(TRUNCATION+1) latitudes
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_BASIC() RESULT(RET) BIND(C)

  CALL SETUP_TEST

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_BASIC

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with two resolutions (regular grid)
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_MULTIPLE_RESOLUTIONS() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION_1 = 39, TRUNCATION_2 = 79

  CALL SETUP_TEST(LCALL_SETUP_TRANS0=.FALSE.)
  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC, KMAX_RESOL=2)

  CALL SETUP_TRANS(KSMAX=TRUNCATION_1, KDGL=2 * (TRUNCATION_1 + 1))
  CALL SETUP_TRANS(KSMAX=TRUNCATION_2, KDGL=2 * (TRUNCATION_2 + 1))

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_MULTIPLE_RESOLUTIONS

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with an odd number of latitudes - should fail
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_ODD_NDGL() RESULT(RET) BIND(C)

  CALL SETUP_TEST

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL - 1)

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_ODD_NDGL

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with octahedral grid of 2*(TRUNCATION+1) latitudes 
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_OCTAHEDRAL() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I

  CALL SETUP_TEST

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_OCTAHEDRAL

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with LDSPLIT option enabled
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_LDSPLIT() RESULT(RET) BIND(C)

  CALL SETUP_TEST

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, LDSPLIT=.TRUE.)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_LDSPLIT

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with stretch factor passed
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_STRETCHING() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I

  CALL SETUP_TEST

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, PSTRET=2.0_JPRD)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_STRETCHING

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with fast Legendre transform 
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_FLT() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I

  CALL SETUP_TEST

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, LDUSEFLT=.TRUE.)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_FLT

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with all fields passed to FFTW at ocne
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_ALL_FFTW() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I

  CALL SETUP_TEST

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, LD_ALL_FFTW=.TRUE.)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_ALL_FFTW

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with Belusov algorithm
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_BELUSOV() RESULT(RET) BIND(C)
  INTEGER(KIND=JPIM) :: ILOEN(NDGL)
  INTEGER(KIND=JPIM) :: I

  CALL SETUP_TEST

  ! Define octahedral grid
  DO I = 1, TRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(NDGL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, KLOEN=ILOEN, LDUSERPNM=.TRUE.)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_BELUSOV

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with LGRIDONLY option enabled
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_GRIDONLY() RESULT(RET) BIND(C)

  CALL SETUP_TEST

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, LDGRIDONLY=.TRUE.)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_GRIDONLY

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with LDSPSETUPONLY option enabled
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_SPSETUPONLY() RESULT(RET) BIND(C)

  CALL SETUP_TEST

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, LDSPSETUPONLY=.TRUE.)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_SPSETUPONLY

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS with LDPNMONLY option enabled
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_PNMONLY() RESULT(RET) BIND(C)

  CALL SETUP_TEST

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, LDPNMONLY=.TRUE.)

  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_PNMONLY

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS works with Legendre polynomials written to and read from disk
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_WRITE_READ_LEGPOLY() RESULT(RET) BIND(C)
  USE PARKIND1, ONLY: JPRB, JPRD
  USE MPL_MODULE, ONLY: MPL_END
  USE OML_MOD, ONLY: OML_MAX_THREADS

  INTEGER(KIND=JPIM) :: IPROC_TEST_NAME
  CHARACTER(LEN=20) :: CLLEGPOLFNAME

#include "abor1.intfb.h"

  CALL SETUP_TEST(LCALL_SETUP_TRANS0=.FALSE.)
  IF (NPROC > 1) THEN
    CALL ABOR1("ECTRANS_TEST_TRANS_API_SETUP_TRANS_WRITE_READ_LEGPOLY: Test does not support" // &
      &        " running with > 1 tasks")
  ENDIF

  ! Build file name for Legendre polynomials, unique to this test
  IPROC_TEST_NAME = MERGE(1, 0, LUSE_MPI)
  IF (JPRB == JPRD) THEN
    WRITE(CLLEGPOLFNAME, "(A,I1,A,I1)") "legpols_dp_mpi", IPROC_TEST_NAME, "xomp", OML_MAX_THREADS()
  ELSE
    WRITE(CLLEGPOLFNAME, "(A,I1,A,I1)") "legpols_sp_mpi", IPROC_TEST_NAME, "xomp", OML_MAX_THREADS()
  ENDIF

  ! First write the polynomial to disk
  CALL SETUP_TRANS0
  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, CDIO_LEGPOL="WRITEF", CDLEGPOLFNAME=CLLEGPOLFNAME)
  CALL TRANS_END

  ! Now setup, reading Legendre polynomials from disk
  CALL SETUP_TRANS0(KPRGPNS=NPROC)
  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, CDIO_LEGPOL="READF", CDLEGPOLFNAME=CLLEGPOLFNAME)

  CALL TRANS_END

  IF (LUSE_MPI) THEN
    CALL MPL_END(LDMEMINFO=.FALSE.)
  ENDIF

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_SETUP_TRANS_WRITE_READ_LEGPOLY

!---------------------------------------------------------------------------------------------------

END MODULE SETUP_TRANS_TEST_SUITE
