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
INTEGER FUNCTION API_TEST_SETUP_TRANS0_EQ_REGIONS() RESULT(RET) BIND(C)
  USE UTIL, ONLY: DETECT_MPIRUN
  USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_END
  USE EC_PARKIND, ONLY: JPIM

  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: NPROC

  LUSE_MPI = DETECT_MPIRUN()

  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    NPROC = MPL_NPROC()
  ELSE
    NPROC = 1
  ENDIF

  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, LDEQ_REGIONS=.TRUE., KPRGPNS=NPROC, KPRGPEW=1, &
    &               KPRTRW=NPROC)

  IF (LUSE_MPI) THEN
    CALL MPL_END(LDMEMINFO=.FALSE.)
  ENDIF

  RET = 0
END FUNCTION API_TEST_SETUP_TRANS0_EQ_REGIONS

!---------------------------------------------------------------------------------------------------

END MODULE SETUP_TRANS0_TEST_SUITE
