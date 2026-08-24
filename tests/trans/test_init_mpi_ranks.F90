! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Test that SETUP_TRANS0 and SETUP_TRANS can be called twice, first on the full MPI_COMM_WORLD and
! then on a reduced communicator holding only the first (NPROC - 1) ranks, exactly as done in
! redistribute.F90 (in ifs-source).
!
! All ranks participate in both MPL_LOCOMM_CREATE calls (collective). Only the subset of ranks that
! belong to the reduced communicator enter the second SETUP_TRANS0 / SETUP_TRANS block.
!
! After the second SETUP_TRANS, an INV_TRANS call exercises TRMTOL which calls
! MPL_ALLTOALLV(KCOMM=MPL_ALL_MS_COMM).  Before FIAT was fixed
! (https://github.com/ecmwf-ifs/fiat/pull/118) this would abort because LGROUPSETUP=.TRUE. prevents
! MPL_GROUPS_CREATE from recreating the communicators for the new NPRTRW, so MPL_ALL_MS_COMM still
! holds Phase-1's larger set of tasks while KRECVCOUNTS is sized for Phase-2's smaller NPRTRW.

PROGRAM TEST_INIT_MPI_RANKS

USE PARKIND1,        ONLY : JPRB, JPIM
USE MPL_MODULE,      ONLY : MPL_INIT, MPL_END, MPL_MYRANK, MPL_NPROC, &
                           & MPL_LOCOMM_CREATE, MPL_SETDFLT_COMM, &
                           & MPL_COMM_OML, MPL_NUMPROC
USE MPL_MPIF,        ONLY : MPI_COMM_WORLD
USE ABORT_TRANS_MOD, ONLY : ABORT_TRANS

IMPLICIT NONE

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "inv_trans.h"
#include "trans_end.h"

INTEGER(KIND=JPIM), PARAMETER :: TRUNCATION = 79
INTEGER(KIND=JPIM), PARAMETER :: NDGL       = 2*(TRUNCATION + 1)

INTEGER(KIND=JPIM) :: NPROC_INIT, MYPROC
INTEGER(KIND=JPIM) :: ICOMM, ICOMM_OLD
INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT
INTEGER(KIND=JPIM) :: NPRGPNS, NPRGPEW, NPRTRW
INTEGER(KIND=JPIM) :: NPROC_REDUCED

! Phase-2 INV_TRANS trigger variables
INTEGER(KIND=JPIM) :: MYSETV2, NLOC_SP, NPROMA2, NGPBLKS2
INTEGER(KIND=JPIM) :: IVSET1(1)
REAL(KIND=JPRB),    ALLOCATABLE :: ZSP2D(:,:), ZGP2D(:,:,:)

CALL MPL_INIT()
NPROC_INIT = MPL_NPROC()
MYPROC     = MPL_MYRANK()

IF (NPROC_INIT < 2) THEN
  CALL ABORT_TRANS("TEST_INIT_MPI_RANKS requires at least 2 MPI tasks")
END IF

! ==================================================================
! Phase 1 – all NPROC_INIT ranks
! ==================================================================
IF (MYPROC == 1) PRINT *, "=== Phase 1: NPROC =", NPROC_INIT

CALL COMPUTE_DECOMP(NPROC_INIT, NPRGPNS, NPRGPEW)
NPRTRW = NPROC_INIT / MAX(NPRGPEW, 1)

CALL SETUP_TRANS0(KPRINTLEV=0, KMAX_RESOL=1, LDMPOFF=.FALSE., &
                  KPRGPNS=NPRGPNS, KPRGPEW=NPRGPEW, KPRTRW=NPRTRW, &
                  LDEQ_REGIONS=.TRUE.)

CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, LDSPLIT=.TRUE.)

CALL TRANS_INQ(KSPEC2=NSPEC2, KGPTOT=NGPTOT)

IF (NSPEC2 <= 0 .OR. NGPTOT <= 0) &
  CALL ABORT_TRANS("Phase 1: TRANS_INQ returned invalid sizes")

PRINT *, "Phase 1 rank", MYPROC, ": KSPEC2 =", NSPEC2, ", KGPTOT =", NGPTOT

CALL TRANS_END()

! ==================================================================================================
! Phase 2 – first (NPROC_INIT - 1) ranks only
!   MPL_LOCOMM_CREATE is collective: ALL ranks call it. Only ranks 1 .. NPROC_REDUCED activate the
!   new communicator and perform the ecTrans setup, exactly as in redistribute.F90 in ifs-source.
! ==================================================================================================
NPROC_REDUCED = NPROC_INIT - 1

CALL MPL_LOCOMM_CREATE(NPROC_REDUCED, ICOMM)

IF (MYPROC <= NPROC_REDUCED) THEN

  CALL MPL_SETDFLT_COMM(ICOMM, ICOMM_OLD)

  IF (MYPROC == 1) PRINT *, "=== Phase 2: NPROC =", NPROC_REDUCED

  CALL COMPUTE_DECOMP(NPROC_REDUCED, NPRGPNS, NPRGPEW)
  NPRTRW = NPROC_REDUCED / MAX(NPRGPEW, 1)

  CALL SETUP_TRANS0(KPRINTLEV=0, KMAX_RESOL=1, LDMPOFF=.FALSE., &
                    KPRGPNS=NPRGPNS, KPRGPEW=NPRGPEW, KPRTRW=NPRTRW, &
                    LDEQ_REGIONS=.TRUE.)

  CALL SETUP_TRANS(KSMAX=TRUNCATION, KDGL=NDGL, LDSPLIT=.TRUE.)

  CALL TRANS_INQ(KSPEC2=NSPEC2, KGPTOT=NGPTOT, KMYSETV=MYSETV2)

  IF (NSPEC2 <= 0 .OR. NGPTOT <= 0) &
    CALL ABORT_TRANS("Phase 2: TRANS_INQ returned invalid sizes")

  PRINT *, "Phase 2 rank", MYPROC, ": KSPEC2 =", NSPEC2, ", KGPTOT =", NGPTOT

  ! Trigger TRMTOL -> MPL_ALLTOALLV(KCOMM=MPL_ALL_MS_COMM).
  ! Without the fiat fix this aborts: SIZE(KRECVCOUNTS)=NPRTRW_phase2 is
  ! smaller than the number of tasks in MPL_ALL_MS_COMM which was built
  ! for NPRTRW_phase1 and was never rebuilt because LGROUPSETUP=.TRUE.
  IVSET1(1) = 1                          ! single field lives on b-set 1
  NLOC_SP   = MERGE(1, 0, MYSETV2 == 1) ! local spectral count on this rank
  NPROMA2   = NGPTOT                     ! one block per task (no KPROMA blocking)
  NGPBLKS2  = 1
  ALLOCATE(ZSP2D(NLOC_SP, NSPEC2))
  ALLOCATE(ZGP2D(NPROMA2, 1, NGPBLKS2))
  ZSP2D = 0.0_JPRB
  ZGP2D = 0.0_JPRB
  IF (MYPROC == 1) PRINT *, "=== Phase 2: calling INV_TRANS (triggers TRMTOL -> MPL_ALL_MS_COMM)"
  CALL INV_TRANS(PSPSC2=ZSP2D, PGP2=ZGP2D, KPROMA=NPROMA2, KVSETSC2=IVSET1)
  DEALLOCATE(ZSP2D, ZGP2D)

  CALL TRANS_END()

END IF

! Restore global communicator before shutdown (as in redistribute.F90)
MPL_COMM_OML(1) = MPI_COMM_WORLD
MPL_NUMPROC     = NPROC_INIT

CALL MPL_END()

IF (MYPROC == 1) PRINT *, "TEST_INIT_MPI_RANKS: PASSED"

CONTAINS

! Compute a 2-D grid decomposition: NPRGPNS x NPRGPEW = N,
! with NPRGPNS >= NPRGPEW, matching the pattern in redistribute.F90.
SUBROUTINE COMPUTE_DECOMP(N, PRGPNS, PRGPEW)
  IMPLICIT NONE
  INTEGER(KIND=JPIM), INTENT(IN)  :: N
  INTEGER(KIND=JPIM), INTENT(OUT) :: PRGPNS, PRGPEW
  INTEGER(KIND=JPIM) :: ISQR, JA, IB

  PRGPNS = N
  PRGPEW = 1
  ISQR = INT(SQRT(REAL(N, JPRB)))
  DO JA = ISQR, N
    IB = N / JA
    IF (JA * IB == N) THEN
      PRGPNS = MAX(JA, IB)
      PRGPEW = MIN(JA, IB)
      EXIT
    END IF
  END DO
END SUBROUTINE COMPUTE_DECOMP

END PROGRAM TEST_INIT_MPI_RANKS
