! (C) Copyright 2026- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM TEST_FLT_BELUSOV

USE PARKIND1, ONLY: JPIM
USE MPL_MODULE, ONLY: MPL_INIT, MPL_END

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_end.h"

INTEGER(KIND=JPIM), PARAMETER :: NT = 383, NL = NT+1

! A linear regular Gaussian grid exposes the small polar Belusov coefficients
! that previously made FLT's interpolation QR return rank zero and abort.

CALL MPL_INIT
CALL SETUP_TRANS0(LDMPOFF=.TRUE.)
CALL SETUP_TRANS(KSMAX=NT,KDGL=NL,LDUSEFLT=.TRUE.,LDUSERPNM=.TRUE.,LDPNMONLY=.TRUE.)
CALL TRANS_END
CALL MPL_END(LDMEMINFO=.FALSE.)

END PROGRAM TEST_FLT_BELUSOV
