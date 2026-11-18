! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SPECNORM_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END
USE UTIL, ONLY: DETECT_MPIRUN

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "specnorm.h"
#include "trans_end.h"

INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 79
INTEGER(KIND=JPIM), PARAMETER :: JPNDGL = 2 * (JPTRUNCATION + 1)
REAL(KIND=JPRB), PARAMETER :: TOLERANCE = 100.0_JPRB * EPSILON(1.0_JPRB)

CONTAINS

!---------------------------------------------------------------------------------------------------

! Test SPECNORM by initialising for each task a local spectral array with each coefficient
! set to the complex value (3 + 4i), except for the m=0 mode which is (3 + 0i). The resulting norm
! is compared with the expected value, taking into account the contributions from the m=0
! coefficients, as well as the m>0 and m<0 coefficients, the latter of which are implicitly stored
! assuming conjugate symmetry.
INTEGER FUNCTION TEST_SPECNORM(LDSPSETUPONLY) RESULT(RET)
  LOGICAL, INTENT(IN) :: LDSPSETUPONLY

  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: IPROC, IMY_PROC, ISPEC2, INUMP, JM, IM, IN, ISP
  INTEGER(KIND=JPIM), ALLOCATABLE :: IMYMS(:)
  INTEGER(KIND=JPIM) :: IASM0(0:JPTRUNCATION)
  REAL(KIND=JPRB), ALLOCATABLE :: ZSPEC(:,:)
  REAL(KIND=JPRB) :: ZNORM(2), ZEXPECTED(2)

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
  CALL SETUP_TRANS(KSMAX=JPTRUNCATION, KDGL=JPNDGL, LDSPSETUPONLY=LDSPSETUPONLY)
  CALL TRANS_INQ(KSPEC2=ISPEC2, KNUMP=INUMP, KASM0=IASM0)
  ALLOCATE(IMYMS(INUMP))
  CALL TRANS_INQ(KMYMS=IMYMS)

  ALLOCATE(ZSPEC(2,ISPEC2))
  ZSPEC(:,:) = 0.0_JPRB
  DO JM = 1, INUMP
    IM = IMYMS(JM)
    DO IN = IM, JPTRUNCATION
      ISP = IASM0(IM) + 2 * (IN - IM)
      ZSPEC(1,ISP) = 3.0_JPRB
      IF (IM > 0) ZSPEC(1,ISP+1) = 4.0_JPRB
    ENDDO
  ENDDO
  ZSPEC(2,:) = -2.0_JPRB * ZSPEC(1,:)

  CALL SPECNORM(PNORM=ZNORM, PSPEC=ZSPEC, KMASTER=1)

  ! T+1 real m=0 modes contribute 9 each; T*(T+1)/2 complex m>0 modes contribute 2*(9+16).
  ZEXPECTED(1) = SQRT(REAL((JPTRUNCATION + 1) * (25 * JPTRUNCATION + 9), JPRB))
  ZEXPECTED(2) = 2.0_JPRB * ZEXPECTED(1)
  RET = 0
  IF (IMY_PROC == 1) THEN
    IF (.NOT. ALL(ABS(ZNORM - ZEXPECTED) <= TOLERANCE * ZEXPECTED)) THEN
      WRITE(*,*) 'SPECNORM: LDSPSETUPONLY=', LDSPSETUPONLY, ', expected ', ZEXPECTED, ', got ', ZNORM
      RET = 1
    ENDIF
  ENDIF
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  DEALLOCATE(IMYMS, ZSPEC)
  CALL TRANS_END
  IF (LUSE_MPI) CALL MPL_END(LDMEMINFO=.FALSE.)
END FUNCTION TEST_SPECNORM

!---------------------------------------------------------------------------------------------------

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SPECNORM() RESULT(RET) BIND(C)
  RET = TEST_SPECNORM(LDSPSETUPONLY=.FALSE.)
END FUNCTION ECTRANS_TEST_TRANS_API_SPECNORM

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_SPECNORM_SPSETUPONLY() RESULT(RET) BIND(C)
  RET = TEST_SPECNORM(LDSPSETUPONLY=.TRUE.)
END FUNCTION ECTRANS_TEST_TRANS_API_SPECNORM_SPSETUPONLY

END MODULE SPECNORM_TEST_SUITE
