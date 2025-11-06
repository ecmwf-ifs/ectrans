! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SETUP_TRANS0_TEST_SUITE

IMPLICIT NONE

#include "setup_trans0.h"

CONTAINS

!---------------------------------------------------------------------------------------------------

! NOTE: SETUP_TRANS0 hardly has any error-prone code - it mostly just sets module variables
! Hence, the only logic we test here is the equal regions initialisation
! For the real tests of initialisation, check the suite for SETUP_TRANS

!---------------------------------------------------------------------------------------------------

! Test SETUP_TRANS0 with equal regions enabled
INTEGER FUNCTION UNIT_TEST_SETUP_TRANS0_EQ_REGIONS() RESULT(RET) BIND(C)
    CALL SETUP_TRANS0(LDMPOFF=.TRUE., LDEQ_REGIONS=.TRUE.)
    RET = 0
END FUNCTION UNIT_TEST_SETUP_TRANS0_EQ_REGIONS

!---------------------------------------------------------------------------------------------------

END MODULE SETUP_TRANS0_TEST_SUITE
