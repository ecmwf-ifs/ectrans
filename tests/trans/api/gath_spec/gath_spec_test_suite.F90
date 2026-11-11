! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE GATH_SPEC_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_ALLREDUCE, MPL_END
USE UTIL, ONLY: DETECT_MPIRUN

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "gath_spec.h"
#include "trans_end.h"

INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 79
INTEGER(KIND=JPIM), PARAMETER :: JPNDGL = 2 * (JPTRUNCATION + 1)

CONTAINS

!---------------------------------------------------------------------------------------------------

! Test GATH_SPEC by initialising on each task a spectral array containing two fields: one storing
! the global index of that spectral coefficient and another with that value multiplied by -2. These
! arrays are then gathered to the last task in the world communicator which checks that the values
! are as expected.
INTEGER FUNCTION TEST_GATH_SPEC(LDSPSETUPONLY) RESULT(RET)
  LOGICAL, INTENT(IN) :: LDSPSETUPONLY

  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: IPROC, IMY_PROC, ISPEC2, ISPEC2G, INUMP, JM, IM, IN, ISP, ISPG
  INTEGER(KIND=JPIM), ALLOCATABLE :: IMYMS(:)
  INTEGER(KIND=JPIM) :: IASM0(0:JPTRUNCATION)
  REAL(KIND=JPRB), ALLOCATABLE :: ZSPECG(:,:), ZSPEC(:,:), ZEXPECTED(:,:)

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
  CALL TRANS_INQ(KSPEC2=ISPEC2, KSPEC2G=ISPEC2G, KNUMP=INUMP, KASM0=IASM0)
  ALLOCATE(IMYMS(INUMP))
  CALL TRANS_INQ(KMYMS=IMYMS)

  ALLOCATE(ZSPECG(2,ISPEC2G), ZSPEC(2,ISPEC2), ZEXPECTED(2,ISPEC2G))
  DO JM = 1, INUMP
    IM = IMYMS(JM)
    DO IN = IM, JPTRUNCATION
      ISP = IASM0(IM) + 2 * (IN - IM)
      ISPG = IM * (2 * JPTRUNCATION + 3 - IM) + 2 * (IN - IM) + 1
      ZSPEC(1,ISP:ISP+1) = [REAL(ISPG, JPRB), REAL(ISPG+1, JPRB)]
    ENDDO
  ENDDO
  ZSPEC(2,:) = -2.0_JPRB * ZSPEC(1,:)

  ! Preserve the deliberately nonzero m=0 imaginary entries to check every copied coefficient.
  CALL GATH_SPEC(PSPECG=ZSPECG, KFGATHG=2, KTO=[IPROC, IPROC], PSPEC=ZSPEC, LDZA0IP=.FALSE.)

  RET = 0
  IF (IMY_PROC == IPROC) THEN
    DO ISPG = 1, ISPEC2G
      ZEXPECTED(1,ISPG) = REAL(ISPG, JPRB)
    ENDDO
    ZEXPECTED(2,:) = -2.0_JPRB * ZEXPECTED(1,:)
    IF (.NOT. ALL(ZSPECG == ZEXPECTED)) THEN
      WRITE(*,*) 'GATH_SPEC: incorrect coefficients on rank ', IMY_PROC, ', LDSPSETUPONLY=', LDSPSETUPONLY
      RET = 1
    ENDIF
  ENDIF
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(RET, CDOPER="MAX")

  DEALLOCATE(IMYMS, ZSPECG, ZSPEC, ZEXPECTED)
  CALL TRANS_END
  IF (LUSE_MPI) CALL MPL_END(LDMEMINFO=.FALSE.)
END FUNCTION TEST_GATH_SPEC

!---------------------------------------------------------------------------------------------------

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_GATH_SPEC() RESULT(RET) BIND(C)
  RET = TEST_GATH_SPEC(LDSPSETUPONLY=.FALSE.)
END FUNCTION ECTRANS_TEST_TRANS_API_GATH_SPEC

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_GATH_SPEC_SPSETUPONLY() RESULT(RET) BIND(C)
  RET = TEST_GATH_SPEC(LDSPSETUPONLY=.TRUE.)
END FUNCTION ECTRANS_TEST_TRANS_API_GATH_SPEC_SPSETUPONLY

END MODULE GATH_SPEC_TEST_SUITE
