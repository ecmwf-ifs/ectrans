! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! Test if we can use TRANS_INQ to access properties initialised through ESETUP_TRANS.

PROGRAM TEST_ETRANS_TRANS_INQ

USE PARKIND1, ONLY: JPIM, JPRB
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS

IMPLICIT NONE

#include "setup_trans0.h"
#include "esetup_trans.h"
#include "trans_inq.h"
#include "etrans_end.h"

INTEGER(KIND=JPIM), PARAMETER :: JP_LAT = 128
INTEGER(KIND=JPIM), PARAMETER :: JP_LON = 128
INTEGER(KIND=JPIM), PARAMETER :: JP_SMAX = JP_LAT / 2 - 1
INTEGER(KIND=JPIM), PARAMETER :: JP_MSMAX = JP_LON / 2 - 1

INTEGER(KIND=JPIM) :: IRESOL, IGPTOTG

CALL SETUP_TRANS0(LDMPOFF=.TRUE.)
CALL ESETUP_TRANS(KSMAX=JP_SMAX, KMSMAX=JP_MSMAX, KDGL=JP_LAT, KDGUX=JP_LAT, KLOEN=[JP_LON], &
  &               KRESOL=IRESOL, PEXWN=1.0_JPRB, PEYWN=1.0_JPRB)

CALL TRANS_INQ(KRESOL=IRESOL, KGPTOTG=IGPTOTG)

CALL ETRANS_END

IF (IGPTOTG /= (JP_LAT * JP_LON)) THEN
  CALL ABORT_TRANS("IGPTOTG does not match expected value")
ELSE
  WRITE(*,'(A,I0)') "IGPTOTG matches expected value = ", IGPTOTG
END IF

END PROGRAM TEST_ETRANS_TRANS_INQ
