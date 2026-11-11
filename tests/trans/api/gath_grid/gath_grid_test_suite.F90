! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE GATH_GRID_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END
USE UTIL, ONLY: DETECT_MPIRUN

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "gath_grid.h"
#include "trans_end.h"

INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 15
INTEGER(KIND=JPIM), PARAMETER :: JPNDGL = 2 * (JPTRUNCATION + 1)
INTEGER(KIND=JPIM), PARAMETER :: JPPROMA = 17

CONTAINS

!---------------------------------------------------------------------------------------------------

! Test GATH_GRID by initialising on each task a gridded array containing two fields: one storing the
! global index of that grid point and another with that value multiplied by -2. These arrays are
! then gathered to the last task in the world communicator which checks that the values are as
! expected. Note that the NPROMA tail of each block is deliberately "poisoned" with -999 to make
! sure that the tails are properly stripped out during the gather operation.
INTEGER FUNCTION TEST_GATH_GRID(LDGRIDONLY) RESULT(RET)
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
  ! Padding must not appear in the gathered field.
  ZGP(:,:,:) = -999.0_JPRB
  DO JGP = 1, IGPTOT
    IROW = MOD(JGP - 1, JPPROMA) + 1
    IBL = (JGP - 1) / JPPROMA + 1
    ZGP(IROW,1,IBL) = REAL(IOFF+JGP, JPRB)
    ZGP(IROW,2,IBL) = -2.0_JPRB * REAL(IOFF+JGP, JPRB)
  ENDDO

  CALL GATH_GRID(PGPG=ZGPG, KPROMA=JPPROMA, KFGATHG=2, KTO=[IPROC, IPROC], PGP=ZGP)

  RET = 0
  IF (IMY_PROC == IPROC) THEN
    DO JGP = 1, IGPTOTG
      ZEXPECTED = [REAL(JGP, JPRB), -2.0_JPRB * REAL(JGP, JPRB)]
      IF (.NOT. ALL(ZGPG(JGP,:) == ZEXPECTED)) THEN
        WRITE(*,*) 'GATH_GRID: incorrect point ', JGP, ' on rank ', IMY_PROC, ', LDGRIDONLY=', LDGRIDONLY
        RET = 1
        EXIT
      ENDIF
    ENDDO
  ENDIF
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  DEALLOCATE(IFRSTLAT, ZGPG, ZGP)
  CALL TRANS_END
  IF (LUSE_MPI) CALL MPL_END(LDMEMINFO=.FALSE.)
END FUNCTION TEST_GATH_GRID

!---------------------------------------------------------------------------------------------------

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_GATH_GRID() RESULT(RET) BIND(C)
  RET = TEST_GATH_GRID(LDGRIDONLY=.FALSE.)
END FUNCTION ECTRANS_TEST_TRANS_API_GATH_GRID

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_GATH_GRID_GRIDONLY() RESULT(RET) BIND(C)
  RET = TEST_GATH_GRID(LDGRIDONLY=.TRUE.)
END FUNCTION ECTRANS_TEST_TRANS_API_GATH_GRID_GRIDONLY

END MODULE GATH_GRID_TEST_SUITE
