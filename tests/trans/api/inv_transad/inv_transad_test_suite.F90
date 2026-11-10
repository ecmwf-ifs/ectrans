! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE INV_TRANSAD_TEST_SUITE

USE PARKIND1, ONLY: JPIM, JPRB, JPRD
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_ALLREDUCE, MPL_END
USE UTIL, ONLY: DETECT_MPIRUN

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "inv_trans.h"
#include "inv_transad.h"
#include "trans_end.h"

! An odd, non-power-of-two block length exercises multiple blocks and a partially
! filled final block, helping expose indexing or alignment assumptions. Dot
! products exclude final-block padding. This tests the blocked API layout;
! the GPU spectral-update race is exposed by the truncations, not by JPPROMA.
INTEGER(KIND=JPIM), PARAMETER :: JPPROMA = 17
INTEGER(KIND=JPIM), PARAMETER :: JPNFIELDS = 2

CONTAINS

!---------------------------------------------------------------------------------------------------

! Test the INV_TRANS/INV_TRANSAD gradient adjoint identity for several truncations. Only the
! north-south and east-west derivative outputs are included in the inner product; scalar values
! are zeroed first so their larger magnitude cannot mask derivative errors.
INTEGER FUNCTION TEST_GRAD_ADJOINT(KTRUNCATION) RESULT(RET)
  INTEGER(KIND=JPIM), INTENT(IN) :: KTRUNCATION

  LOGICAL :: LUSE_MPI
  INTEGER(KIND=JPIM) :: IPROC, ISPEC2, INUMP, IGPTOT, IGPBLKS, JM, IM, IN, ISP, JGP, IROW, IBL
  INTEGER(KIND=JPIM), ALLOCATABLE :: IMYMS(:), IASM0(:)
  REAL(KIND=JPRB), ALLOCATABLE :: ZSPEC(:,:), ZADJ(:,:), ZGP(:,:,:)
  REAL(KIND=JPRD) :: ZDOT(2), ZWEIGHT, ZERROR

  LUSE_MPI = DETECT_MPIRUN()
  IF (LUSE_MPI) THEN
    CALL MPL_INIT
    IPROC = MPL_NPROC()
  ELSE
    IPROC = 1
  ENDIF

  CALL SETUP_TRANS0(LDMPOFF=.NOT. LUSE_MPI, KPRGPNS=IPROC, KPRGPEW=1, KPRTRW=IPROC)
  CALL SETUP_TRANS(KSMAX=KTRUNCATION, KDGL=2*(KTRUNCATION+1))
  ALLOCATE(IASM0(0:KTRUNCATION))
  CALL TRANS_INQ(KSPEC2=ISPEC2, KNUMP=INUMP, KGPTOT=IGPTOT, KASM0=IASM0)
  ALLOCATE(IMYMS(INUMP))
  CALL TRANS_INQ(KMYMS=IMYMS)

  IGPBLKS = (IGPTOT - 1) / JPPROMA + 1
  ALLOCATE(ZSPEC(JPNFIELDS,ISPEC2), ZADJ(JPNFIELDS,ISPEC2), ZGP(JPPROMA,3*JPNFIELDS,IGPBLKS))
  ZSPEC(:,:) = 0.0_JPRB
  DO JM = 1, INUMP
    IM = IMYMS(JM)
    DO IN = IM, KTRUNCATION
      ISP = IASM0(IM) + 2*(IN-IM)
      ZSPEC(1,ISP) = REAL(IM+1, JPRB) / REAL(IN+1, JPRB)
      IF (IM > 0) ZSPEC(1,ISP+1) = -0.5_JPRB * ZSPEC(1,ISP)
    ENDDO
  ENDDO
  ZSPEC(2,:) = -2.0_JPRB * ZSPEC(1,:)
  ZGP(:,:,:) = 0.0_JPRB
  CALL INV_TRANS(PSPSCALAR=ZSPEC, PGP=ZGP, KPROMA=JPPROMA, LDSCDERS=.TRUE.)

  ! Exclude scalar values so their much larger magnitude cannot hide derivative errors.
  ZGP(:,1:JPNFIELDS,:) = 0.0_JPRB
  ZDOT(:) = 0.0_JPRD
  DO JGP = 1, IGPTOT
    IROW = MOD(JGP-1, JPPROMA) + 1
    IBL = (JGP-1) / JPPROMA + 1
    ZDOT(1) = ZDOT(1) + SUM(REAL(ZGP(IROW,:,IBL), JPRD)**2)
  ENDDO

  ZADJ(:,:) = 0.0_JPRB
  CALL INV_TRANSAD(PSPSCALAR=ZADJ, PGP=ZGP, KPROMA=JPPROMA, LDSCDERS=.TRUE.)
  DO JM = 1, INUMP
    IM = IMYMS(JM)
    ZWEIGHT = REAL(1 + MIN(1, IM), JPRD)
    DO IN = IM, KTRUNCATION
      ISP = IASM0(IM) + 2*(IN-IM)
      ZDOT(2) = ZDOT(2) + ZWEIGHT * SUM(REAL(ZSPEC(:,ISP), JPRD) * REAL(ZADJ(:,ISP), JPRD))
      IF (IM > 0) THEN
        ZDOT(2) = ZDOT(2) + ZWEIGHT * SUM(REAL(ZSPEC(:,ISP+1), JPRD) * REAL(ZADJ(:,ISP+1), JPRD))
      ENDIF
    ENDDO
  ENDDO
  IF (LUSE_MPI) CALL MPL_ALLREDUCE(ZDOT, CDOPER="SUM", LDREPROD=.FALSE.)

  RET = 0
  ZERROR = ABS(ZDOT(1) - ZDOT(2)) / MAX(ABS(ZDOT(1)), TINY(1.0_JPRD))
  IF (.NOT. (ZDOT(1) > 0.0_JPRD .AND. ZERROR <= 1000.0_JPRD * EPSILON(1.0_JPRB))) THEN
    WRITE(*,*) 'INV_TRANS gradient adjoint: truncation=', KTRUNCATION, ', dot products=', ZDOT, &
      & ', relative error=', ZERROR
    RET = 1
  ENDIF

  DEALLOCATE(IMYMS, IASM0, ZSPEC, ZADJ, ZGP)
  CALL TRANS_END
  IF (LUSE_MPI) CALL MPL_END(LDMEMINFO=.FALSE.)
END FUNCTION TEST_GRAD_ADJOINT

!---------------------------------------------------------------------------------------------------

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_INV_TRANSAD_GRAD_ADJOINT_T7() RESULT(RET) BIND(C)
  RET = TEST_GRAD_ADJOINT(7)
END FUNCTION ECTRANS_TEST_TRANS_API_INV_TRANSAD_GRAD_ADJOINT_T7

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_INV_TRANSAD_GRAD_ADJOINT_T11() RESULT(RET) BIND(C)
  RET = TEST_GRAD_ADJOINT(11)
END FUNCTION ECTRANS_TEST_TRANS_API_INV_TRANSAD_GRAD_ADJOINT_T11

INTEGER FUNCTION ECTRANS_TEST_TRANS_API_INV_TRANSAD_GRAD_ADJOINT_T14() RESULT(RET) BIND(C)
  RET = TEST_GRAD_ADJOINT(14)
END FUNCTION ECTRANS_TEST_TRANS_API_INV_TRANSAD_GRAD_ADJOINT_T14

END MODULE INV_TRANSAD_TEST_SUITE
