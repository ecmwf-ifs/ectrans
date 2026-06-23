! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE GPNORM_TRANS_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "inv_trans.h"
#include "gpnorm_trans.h"
#include "trans_end.h"

! Spectral truncation used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 79

! Number of latitudes used for all tests below
INTEGER(KIND=JPIM), PARAMETER :: JPNGDL = 2 * (JPTRUNCATION + 1)

! Tolerance for "close to zero"
REAL(KIND=JPRB), PARAMETER :: PPTOLERANCE = 100.0_JPRB * EPSILON(1.0_JPRB)

! JPPROMA blocking factor
INTEGER(KIND=JPIM), PARAMETER :: JPPROMA = 16

LOGICAL :: LUSE_MPI

CONTAINS

!---------------------------------------------------------------------------------------------------

! Approximate equality check for reals
ELEMENTAL LOGICAL FUNCTION APPROX_EQ(PA, PB, PTOL) RESULT(KRET)
  REAL(KIND=JPRB), INTENT(IN) :: PA
  REAL(KIND=JPRB), INTENT(IN) :: PB
  REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: PTOL

  IF (PRESENT(PTOL)) THEN
    KRET = ABS(PA - PB) < PTOL
  ELSE
    KRET = ABS(PA - PB) < PPTOLERANCE
  END IF
END FUNCTION APPROX_EQ

! Setup fixture
SUBROUTINE SETUP_TEST(KSPEC2, KGPBLKS, KMY_PROC)
  USE UTIL, ONLY: DETECT_MPIRUN

  INTEGER(KIND=JPIM), INTENT(OUT) :: KSPEC2
  INTEGER(KIND=JPIM), INTENT(OUT) :: KGPBLKS
  INTEGER(KIND=JPIM), INTENT(OUT) :: KMY_PROC

  INTEGER(KIND=JPIM) :: ILOEN(JPNGDL)
  INTEGER(KIND=JPIM) :: I
  INTEGER(KIND=JPIM) :: IPROC
  INTEGER(KIND=JPIM) :: KGPTOT

  ! Set up MPI
  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    IPROC = MPL_NPROC()
    KMY_PROC = MPL_MYRANK()
  ELSE
    IPROC = 1
    KMY_PROC = 1
  ENDIF

  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=IPROC, KPRGPEW=1, KPRTRW=IPROC)

  ! Define octahedral grid
  DO I = 1, JPTRUNCATION + 1
    ILOEN(I) = 20 + 4 * I
    ILOEN(JPNGDL - I + 1) = ILOEN(I)
  END DO

  CALL SETUP_TRANS(KSMAX=JPTRUNCATION, KDGL=JPNGDL, KLOEN=ILOEN)

  CALL TRANS_INQ(KSPEC2=KSPEC2, KGPTOT=KGPTOT)

  ! Number of JPPROMA blocks
  KGPBLKS = (KGPTOT - 1) / JPPROMA + 1
END SUBROUTINE SETUP_TEST

! Tear down fixture
SUBROUTINE CLEANUP_TEST

  CALL TRANS_END

  IF (LUSE_MPI) THEN
    CALL MPL_END(LDMEMINFO=.FALSE.)
  ENDIF
END SUBROUTINE CLEANUP_TEST

!---------------------------------------------------------------------------------------------------

! Test GPNORM_TRANS average
INTEGER FUNCTION ECTRANS_TEST_TRANS_API_GPNORM_TRANS_AVE() RESULT(KRET) BIND(C)
  INTEGER(KIND=JPIM) :: ISPEC2, IGPBLKS, IMY_PROC
  REAL(KIND=JPRB), ALLOCATABLE :: ZSPECTRAL_FIELD(:,:), ZGRID_POINT_FIELD(:,:,:)
  INTEGER(KIND=JPIM) :: INUM_MY_ZON_WNS
  INTEGER(KIND=JPIM), ALLOCATABLE :: IMY_ZON_WNS(:)
  REAL(KIND=JPRB) :: ZAVE(1), ZMIN(1), ZMAX(1)

  REAL(KIND=JPRB), PARAMETER :: ZAVERAGE = 3.1415926535_JPRB

  CALL SETUP_TEST(ISPEC2, IGPBLKS, IMY_PROC)

  ALLOCATE(ZSPECTRAL_FIELD(1,ISPEC2))
  ALLOCATE(ZGRID_POINT_FIELD(JPPROMA,1,IGPBLKS))
  ZSPECTRAL_FIELD(:,:) = 0.0

  ! Initialise our spectral field array so that the (0, 0) mode has our average, and all other values
  ! are zero
  CALL TRANS_INQ(KNUMP=INUM_MY_ZON_WNS)
  ALLOCATE(IMY_ZON_WNS(INUM_MY_ZON_WNS))
  CALL TRANS_INQ(KMYMS=IMY_ZON_WNS)
  IF (ANY(IMY_ZON_WNS == 0)) THEN
      ZSPECTRAL_FIELD(1,1) = ZAVERAGE
  ENDIF

  ! Transform to grid point space
  CALL INV_TRANS(PSPSCALAR=ZSPECTRAL_FIELD, PGP=ZGRID_POINT_FIELD, KPROMA=JPPROMA)

  ! Now compute the average and check it's correct
  CALL GPNORM_TRANS(PGP=ZGRID_POINT_FIELD, KFIELDS=1, KPROMA=JPPROMA, PAVE=ZAVE, PMIN=ZMIN, &
    &               PMAX=ZMAX, LDAVE_ONLY=.TRUE.)

  IF (IMY_PROC == 1) THEN
    IF (APPROX_EQ(ZAVE(1), ZAVERAGE)) THEN
        KRET = 0
    ELSE
        KRET = 1
    ENDIF
  ELSE
    KRET = 0
  ENDIF

  IF (LUSE_MPI) CALL MPL_ALLREDUCE(KRET, CDOPER="MAX")

  ! Tear down everything
  DEALLOCATE(IMY_ZON_WNS, ZGRID_POINT_FIELD, ZSPECTRAL_FIELD)
  CALL CLEANUP_TEST
END FUNCTION ECTRANS_TEST_TRANS_API_GPNORM_TRANS_AVE

!---------------------------------------------------------------------------------------------------

END MODULE GPNORM_TRANS_TEST_SUITE
