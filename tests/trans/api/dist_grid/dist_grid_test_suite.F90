! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE DIST_GRID_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END
USE UTIL, ONLY: DETECT_MPIRUN

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "dist_grid.h"
#include "trans_end.h"

INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 15
INTEGER(KIND=JPIM), PARAMETER :: JPNDGL = 2 * (JPTRUNCATION + 1)
INTEGER(KIND=JPIM), PARAMETER :: JPPROMA = 17

CONTAINS

!---------------------------------------------------------------------------------------------------

! Test DIST_GRID by initialising a global gridded array on the last task in the world communicator,
! with each grid point containing two fields: one storing a real version of its global index, and
! another with that value multiplied by -2. This field is distributed to all other tasks, who then
! check their values against the expected values derived from the indexing arrays.
INTEGER FUNCTION TEST_DIST_GRID(LDGRIDONLY) RESULT(RET)
  LOGICAL, INTENT(IN) :: LDGRIDONLY

  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: IPROC, IMY_PROC, IGPTOT, IGPTOTG, IGPBLKS, IREGION, IOFF
  INTEGER(KIND=JPIM) :: JLAT, JGP, IROW, IBL
  INTEGER(KIND=JPIM) :: ILOEN(JPNDGL)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IFRSTLAT(:)
  REAL(KIND=JPRB), ALLOCATABLE :: ZGPG(:,:), ZGP(:,:,:)
  REAL(KIND=JPRB) :: ZEXPECTED(2)

  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    IPROC = MPL_NPROC()
    IMY_PROC = MPL_MYRANK()
  ELSE
    IPROC = 1
    IMY_PROC = 1
  ENDIF

  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=IPROC, KPRGPEW=1, KPRTRW=IPROC)
  DO JLAT = 1, JPNDGL / 2
    ILOEN(JLAT) = 20 + 4 * JLAT
    ILOEN(JPNDGL - JLAT + 1) = ILOEN(JLAT)
  ENDDO
  CALL SETUP_TRANS(KSMAX=JPTRUNCATION, KDGL=JPNDGL, KLOEN=ILOEN, LDSPLIT=.FALSE., LDGRIDONLY=LDGRIDONLY)
  ALLOCATE(IFRSTLAT(IPROC))
  CALL TRANS_INQ(KGPTOT=IGPTOT, KGPTOTG=IGPTOTG, KFRSTLAT=IFRSTLAT, KMY_REGION_NS=IREGION)
  ! With one east-west region and no split latitudes, each local grid is a contiguous latitude band.
  IOFF = SUM(ILOEN(1:IFRSTLAT(IREGION)-1))
  IGPBLKS = (IGPTOT - 1) / JPPROMA + 1
  ALLOCATE(ZGPG(IGPTOTG,2), ZGP(JPPROMA,2,IGPBLKS))
  ZGPG(:,:) = 0.0_JPRB
  IF (IMY_PROC == IPROC) THEN
    DO JGP = 1, IGPTOTG
      ZGPG(JGP,1) = REAL(JGP, JPRB)
      ZGPG(JGP,2) = -2.0_JPRB * REAL(JGP, JPRB)
    ENDDO
  ENDIF

  CALL DIST_GRID(PGPG=ZGPG, KPROMA=JPPROMA, KFDISTG=2, KFROM=[IPROC, IPROC], PGP=ZGP)

  RET = 0
  DO JGP = 1, IGPTOT
    IROW = MOD(JGP - 1, JPPROMA) + 1
    IBL = (JGP - 1) / JPPROMA + 1
    ZEXPECTED = [REAL(IOFF+JGP, JPRB), -2.0_JPRB * REAL(IOFF+JGP, JPRB)]
    IF (.NOT. ALL(ZGP(IROW,:,IBL) == ZEXPECTED)) THEN
      WRITE(*,*) 'DIST_GRID: incorrect point ', JGP, ' on rank ', IMY_PROC, ', LDGRIDONLY=', LDGRIDONLY
      RET = 1
      EXIT
    ENDIF
  ENDDO
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  DEALLOCATE(IFRSTLAT, ZGPG, ZGP)
  CALL TRANS_END
  IF (LUSE_MPI) CALL MPL_END(LDMEMINFO=.FALSE.)
END FUNCTION TEST_DIST_GRID

!---------------------------------------------------------------------------------------------------

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_DIST_GRID() RESULT(RET) BIND(C)
  RET = TEST_DIST_GRID(LDGRIDONLY=.FALSE.)
END FUNCTION ECTRANS_TEST_TRANS_API_DIST_GRID

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_DIST_GRID_GRIDONLY() RESULT(RET) BIND(C)
  RET = TEST_DIST_GRID(LDGRIDONLY=.TRUE.)
END FUNCTION ECTRANS_TEST_TRANS_API_DIST_GRID_GRIDONLY

END MODULE DIST_GRID_TEST_SUITE
