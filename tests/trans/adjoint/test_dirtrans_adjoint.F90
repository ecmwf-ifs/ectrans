! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! ==================================================================================================
! DIR_TRANS adjoint test
! ==================================================================================================
!
! This program performs a rudimentary check of tangent-linear/adjoint correspondence of 
! DIR_TRANS and DIR_TRANSAD.
!
! The program checks the correspondence of <DIR_TRANS(X1), Y2> and <X1, DIR_TRANSAD(Y2)>
! which with infinite precision should match exactly. In practice there is some divergence due to
! rounding errors. In this program we check whether the two expressions are the same to within a
! tolerance of 20000 * machine epsilon.
!
! In this test X1, the "global state vector", is comprised of scalar, and u-v fields
! defined on 9 model levels. The correspondence is computed across the whole state vector.
!
! ==================================================================================================

PROGRAM TEST_DIRTRANS_ADJOINT

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE, ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC, MPL_BARRIER, MPL_END
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
USE UTILS, ONLY: DETECT_MPIRUN, SCALPRODSP, SCALPRODGP

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 159 ! T159 spectral resolution
INTEGER(KIND=JPIM), PARAMETER :: JPPROMA = 16
INTEGER(KIND=JPIM), PARAMETER :: JP_NUMLEVELS_G = 9
INTEGER(KIND=JPIM), PARAMETER :: JPNLAT = 2 * (JPTRUNCATION + 1)

INTEGER(KIND=JPIM) :: INPROC, IMYPROC, IPRGPNS, IPRGPEW, IPRTRW, IPRTRV, IGPTOTG, IGPTOT, IGPBLKS
INTEGER(KIND=JPIM) :: ISPEC2G, ISPEC2
INTEGER(KIND=JPIM) :: INUM_LEVELS
INTEGER(KIND=JPIM) :: IMYSETV
INTEGER(KIND=JPIM) :: IVSET(JP_NUMLEVELS_G)
INTEGER(KIND=JPIM) :: NLOEN(JPNLAT)
INTEGER(KIND=JPIM) :: ITOSP(JP_NUMLEVELS_G), ITOGP(3*JP_NUMLEVELS_G)
INTEGER(KIND=JPIM) :: JLEV, JA, JB

LOGICAL :: LLUSE_MPI
INTEGER(KIND=JPIM) :: IOUT = 6, IERR = 0 ! STDOUT and STDERR
REAL(KIND=JPRB) , ALLOCATABLE :: ZSPECX(:,:), ZSPECY(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZVORX(:,:), ZVORY(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZDIVX(:,:), ZDIVY(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZGX(:,:,:), ZGY(:,:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZSPECXG(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZVORXG(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZDIVXG(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZGXG(:,:)
REAL(KIND=JPRB) :: ADJ_VALUE_1
REAL(KIND=JPRB) :: ADJ_VALUE_2
REAL(KIND=JPRB) :: ZRELATIVE_ERROR
INTEGER(KIND=JPIM) :: N
INTEGER(KIND=JPIM), ALLOCATABLE :: SEED(:)

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "dir_trans.h"
#include "dir_transad.h"
#include "dist_grid.h"
#include "dist_spec.h"
#include "trans_end.h"

! Fix random number seed
CALL RANDOM_SEED(SIZE=N)
ALLOCATE(SEED(N))
SEED(:) = 1
CALL RANDOM_SEED(PUT=SEED)

LLUSE_MPI = DETECT_MPIRUN()

! Set up MPI
IF (LLUSE_MPI) THEN
  CALL MPL_INIT
  IMYPROC = MPL_MYRANK()
  INPROC = MPL_NPROC()
ELSE
  IMYPROC = 1
  INPROC  = 1
ENDIF

! Only output to stdout on first task
IF (INPROC > 1) THEN
  IF (IMYPROC /= 1) THEN
    OPEN(UNIT=IOUT, FILE='/dev/null')
  ENDIF
ENDIF

! Compute E-W and V-W set sizes
DO JA = INT(SQRT(REAL(INPROC,JPRB))), INPROC
  JB = INPROC / JA
  IF (JA * JB == INPROC) THEN
    IPRGPNS = MAX(JA, JB)
    IPRGPEW = MIN(JA, JB)
    IPRTRW  = MAX(JA, JB)
    IPRTRV  = MIN(JA, JB)
  ENDIF
ENDDO

IMYSETV = MOD(IMYPROC - 1, IPRTRV) + 1

! Use a full Gaussian grid
NLOEN(:) = 2*JPNLAT

! Initialise ecTrans
CALL SETUP_TRANS0(KOUT=IOUT, KERR=IERR, KPRGPNS=IPRGPNS, KPRGPEW=IPRGPEW, KPRTRW=IPRTRW, &
  &               LDMPOFF=.NOT. LLUSE_MPI, KPRINTLEV=0)
CALL SETUP_TRANS(KSMAX=JPTRUNCATION, KDGL=JPNLAT, KLOEN=NLOEN, LDSPLIT=.TRUE.)
CALL TRANS_INQ(KSPEC2=ISPEC2, KSPEC2G=ISPEC2G, KGPTOT=IGPTOT, KGPTOTG=IGPTOTG)

IGPBLKS = (IGPTOT-1)/JPPROMA+1

! Determine number of local levels
INUM_LEVELS = 0
DO JLEV = 1, JP_NUMLEVELS_G
  IVSET(JLEV) = MOD(JLEV, IPRTRV) + 1
  IF (IVSET(JLEV) == IMYSETV) THEN
    INUM_LEVELS = INUM_LEVELS + 1
  ENDIF
ENDDO

! Initially task 1 has all the fields
ITOSP(:) = 1
ITOGP(:) = 1

! ===== Allocate and initialize spectral data =====
ALLOCATE(ZSPECXG(JP_NUMLEVELS_G,ISPEC2G))
ALLOCATE(ZVORXG(JP_NUMLEVELS_G,ISPEC2G))
ALLOCATE(ZDIVXG(JP_NUMLEVELS_G,ISPEC2G))

ALLOCATE(ZSPECX(INUM_LEVELS,ISPEC2))
ALLOCATE(ZSPECY(INUM_LEVELS,ISPEC2))
ALLOCATE(ZVORX(INUM_LEVELS,ISPEC2))
ALLOCATE(ZVORY(INUM_LEVELS,ISPEC2))
ALLOCATE(ZDIVX(INUM_LEVELS,ISPEC2))
ALLOCATE(ZDIVY(INUM_LEVELS,ISPEC2))

IF (IMYPROC == 1) THEN
  CALL RANDOM_NUMBER(ZSPECXG)
  ZSPECXG(:,:) = 0.1_JPRB * (1.0_JPRB - 2.0_JPRB * ZSPECXG(:,:))
  CALL RANDOM_NUMBER(ZVORXG)
  ZVORXG(:,:) = 0.1_JPRB * (1.0_JPRB - 2.0_JPRB * ZVORXG(:,:))
  CALL RANDOM_NUMBER(ZDIVXG)
  ZDIVXG(:,:) = 0.1_JPRB * (1.0_JPRB - 2.0_JPRB * ZDIVXG(:,:))
ENDIF

CALL DIST_SPEC(PSPECG=ZSPECXG, KFDISTG=JP_NUMLEVELS_G, KFROM=ITOSP, PSPEC=ZSPECX, KVSET=IVSET)
CALL DIST_SPEC(PSPECG=ZVORXG, KFDISTG=JP_NUMLEVELS_G, KFROM=ITOSP, PSPEC=ZVORX, KVSET=IVSET)
CALL DIST_SPEC(PSPECG=ZDIVXG, KFDISTG=JP_NUMLEVELS_G, KFROM=ITOSP, PSPEC=ZDIVX, KVSET=IVSET)

! ===== Allocate and initialize gridpoint data =====
ALLOCATE(ZGXG(IGPTOTG,3*JP_NUMLEVELS_G))
ALLOCATE(ZGX(JPPROMA,3*JP_NUMLEVELS_G,IGPBLKS))
ALLOCATE(ZGY(JPPROMA,3*JP_NUMLEVELS_G,IGPBLKS))

IF (IMYPROC == 1) THEN
  CALL RANDOM_NUMBER(ZGXG)
  ZGXG(:,:) = (1.0_JPRB-2.0_JPRB*ZGXG(:,:))
ENDIF

CALL DIST_GRID(PGPG=ZGXG, KFDISTG=3*JP_NUMLEVELS_G, KFROM=ITOGP, PGP=ZGX, KPROMA=JPPROMA)

! ===== Compute adjoint dirtrans and gather result on proc 1 =====
! i.e. dirtrans(rgpx) = (rspscalary, rspvory, rspdivy)

CALL DIR_TRANS(PSPSCALAR=ZSPECY, PSPVOR=ZVORY, PSPDIV=ZDIVY, PGP=ZGX, KPROMA=JPPROMA, &
  &            KVSETSC=IVSET, KVSETUV=IVSET)

! ===== Compute: adj_value1 = <dirtrans(rgpx), (rspscalarx, rspvorx, rspdivx)> =====
! i.e. adj_value1 = <(rspscalary, rspvory, rspdivy), (rspscalarx, rspvorx, rspdivx)>

ADJ_VALUE_1 = SCALPRODSP(ZSPECX, ZSPECY, IVSET, INUM_LEVELS, JP_NUMLEVELS_G, ISPEC2, ISPEC2G, JPTRUNCATION, IMYPROC) + &
  & SCALPRODSP(ZVORX, ZVORY, IVSET, INUM_LEVELS, JP_NUMLEVELS_G, ISPEC2, ISPEC2G, JPTRUNCATION, IMYPROC ) + &
  & SCALPRODSP(ZDIVX, ZDIVY, IVSET, INUM_LEVELS, JP_NUMLEVELS_G, ISPEC2, ISPEC2G, JPTRUNCATION, IMYPROC)

! ===== Compute dirtrans_adj and gather result on proc 1 =====
! i.e. dirtrans_adj(rspscalarx, rspvorx, rspdivx) = rgpy

CALL DIR_TRANSAD(PSPSCALAR=ZSPECX, PSPVOR=ZVORX, PSPDIV=ZDIVX, PGP=ZGY, KPROMA=JPPROMA, &
  &              KVSETSC=IVSET, KVSETUV=IVSET)

! ===== Compute: adj_value2 = <rgpx, dirtrans_adj(rspscalarx, rspvorx, rspdivx)> =====
! i.e. adj_value2 = <rgpy, rgpx>

ADJ_VALUE_2 = SCALPRODGP(ZGY, ZGX, JPPROMA, 3 * JP_NUMLEVELS_G, IGPBLKS, IGPTOT, IGPTOTG, IMYPROC)

! Only task 1 should perform the correctness check
IF (IMYPROC == 1) THEN
  ! ===== Compare inner products =====
  ! i.e. <dirtrans(rgpx), (rspscalarx, rspvorx, rspdivx)> == <rgpx, dirtrans_adj(rspscalarx, rspvorx, rspdivx)>

  ZRELATIVE_ERROR = ABS(ADJ_VALUE_1 - ADJ_VALUE_2) / ABS(ADJ_VALUE_1)

  WRITE(IOUT, '(A,1E30.15)') '<Fx,y>  = ', ADJ_VALUE_1
  WRITE(IOUT, '(A,1E30.15)') '<x,F*y> = ', ADJ_VALUE_2
  WRITE(IOUT, '(A,1E20.15)') 'Relative error = ', ZRELATIVE_ERROR

  ! Abort if relative error is > 20000 * machine epsilon
  ! All tested compilers seem to be happy with a threshold of 20000, thought it is a bit arbitrary
  IF (ZRELATIVE_ERROR > 20000.0*EPSILON(1.0_JPRB)) THEN
    WRITE(IOUT, '(A)') '*******************************'
    WRITE(IOUT, '(A)') 'Adjoint test failed'
    WRITE(IOUT, '(A)') 'Relative error greater than 20000 * machine epsilon'
    WRITE(IOUT, '(1E20.15,A3,1E20.15)') ZRELATIVE_ERROR, ' > ', 20000.0*EPSILON(1.0_JPRB)
    WRITE(IOUT, '(A)') '*******************************'
    FLUSH(IOUT)
    CALL TRANS_END
    CALL ABORT_TRANS("Adjoint test failed")
  ENDIF
ENDIF

CALL TRANS_END

IF (LLUSE_MPI) THEN
  CALL MPL_BARRIER()
  CALL MPL_END
ENDIF

END PROGRAM TEST_DIRTRANS_ADJOINT
