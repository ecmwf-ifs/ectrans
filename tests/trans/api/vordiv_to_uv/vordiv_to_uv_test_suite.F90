! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE VORDIV_TO_UV_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB, JPRD

IMPLICIT NONE

#include "setup_trans0.h"
#include "vordiv_to_uv.h"
#include "trans_end.h"

LOGICAL :: LUSE_MPI
INTEGER(KIND=JPIM) :: NPROC

CONTAINS

!---------------------------------------------------------------------------------------------------


! Setup fixture
SUBROUTINE SETUP_TEST()
  USE UTIL, ONLY: DETECT_MPIRUN, ENABLE_FPE
  USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC

  ! Set up MPI
  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    NPROC = MPL_NPROC()
  ELSE
    NPROC = 1
  ENDIF

  CALL ENABLE_FPE() ! Can be disabled by setting ECTRANS_TEST_ENABLE_FPE environment variable to "0"

  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=NPROC, KPRGPEW=1, KPRTRW=NPROC)
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

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_VORDIV_TO_UV_T1() RESULT(RET) BIND(C)

  INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION = 1
  INTEGER(KIND=JPIM), PARAMETER :: NFLD = 1
  INTEGER(KIND=JPIM), PARAMETER :: NSPEC2 = (TRUNCATION + 1) * (TRUNCATION + 2)

  REAL(KIND=JPRB) :: RSPVOR(NFLD, NSPEC2)
  REAL(KIND=JPRB) :: RSPDIV(NFLD, NSPEC2)
  REAL(KIND=JPRB) :: RSPU(NFLD, NSPEC2)
  REAL(KIND=JPRB) :: RSPV(NFLD, NSPEC2)

  RSPVOR = 1.0_JPRB
  RSPDIV = 1.0_JPRB

  CALL SETUP_TEST
  CALL VORDIV_TO_UV(PSPVOR=RSPVOR, PSPDIV=RSPDIV, PSPU=RSPU, PSPV=RSPV, KSMAX=TRUNCATION)

  ! RSPU and RSPV contain NaNs if FPE trapping is not working.
  WRITE(0,*) 'PSPU =', RSPU
  WRITE(0,*) 'PSPV =', RSPV
  CALL END_TEST

  RET = 0
END FUNCTION ECTRANS_TEST_TRANS_API_VORDIV_TO_UV_T1

!---------------------------------------------------------------------------------------------------

END MODULE VORDIV_TO_UV_TEST_SUITE
