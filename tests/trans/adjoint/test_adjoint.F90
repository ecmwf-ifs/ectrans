! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ==================================================================================================
! Adjoint test
! ==================================================================================================
!
! This program performs a rudimentary check of tangent-linear/adjoint correspondence of the inverse
! and direct spectral transform.
!
! The program checks the correspondence of <DIR_TRANS(INV_TRANS(X)), Y> and
! <X, INV_TRANSAD(DIR_TRANSAD(Y))>, which with infinite precision should match exactly. In practice
! there is some divergence due to rounding errors. In this program we check whether the two
! expressions are the same to within a tolerance of 2000 * machine epsilon.
!
! The check is only performed for scalar fields (PSPSCALAR). Wind fields are not checked.
!
! ==================================================================================================

PROGRAM TEST_ADJOINT

USE PARKIND1,        ONLY: JPIM, JPRB
USE MPL_MODULE,      ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC, MPL_BARRIER, MPL_END
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
USE UTILS,           ONLY: DETECT_MPIRUN, SCALPRODSP

IMPLICIT NONE

INTEGER(KIND=JPIM) :: NSMAX, NDGL, NPROC, NPRGPNS, NPRGPEW, NPRTRW, NPRTRV
INTEGER(KIND=JPIM) :: NOUT, NERR, MYPROC, NSPECG, NSPEC2G
INTEGER(KIND=JPIM) :: NFLEV, NFLEVG
INTEGER(KIND=JPIM) :: NSPEC2, NGPTOT, NPROMA, NGPBLKS, MYSETV
INTEGER(KIND=JPIM), ALLOCATABLE :: IVSET(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: NLOEN(:)
INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) , ALLOCATABLE :: ZSPECX(:,:), ZSPECY(:,:), ZSPECP(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZGX(:,:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZSPECYG(:,:), ZSPECXG(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZRANDSP(:)
REAL(KIND=JPRB) :: ZSC1, ZSC2, ZRELATIVE_ERROR
INTEGER(KIND=JPIM) :: JA, JB, I
LOGICAL :: LUSE_MPI
INTEGER(KIND=JPIM) :: N
INTEGER(KIND=JPIM), ALLOCATABLE :: SEED(:)

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "dir_trans.h"
#include "inv_trans.h"
#include "dir_transad.h"
#include "inv_transad.h"
#include "dist_grid.h"
#include "dist_spec.h"
#include "trans_end.h"

! Fix random number seed
CALL RANDOM_SEED(SIZE=N)
ALLOCATE(SEED(N))
SEED(:) = 1
CALL RANDOM_SEED(PUT=SEED)

LUSE_MPI = DETECT_MPIRUN()

NDGL = 32 ! Number of latitudes from pole to equator
NFLEVG = 9 ! Number of levels
NPROMA = 8 ! Gridpoint block size

! Determine spectral space parameters
NSMAX = (2 * NDGL - 1) / 3 ! Full Gaussian grid
NSPECG = (NSMAX+1)*(NSMAX+2)/2
NSPEC2G = NSPECG*2

IF (LUSE_MPI) THEN
  CALL MPL_INIT
  MYPROC = MPL_MYRANK()
  NPROC = MPL_NPROC()
ELSE
  MYPROC = 1
  NPROC  = 1
ENDIF

! STDOUT and STDERR
NOUT = 6
NERR = 0

! Only output to stdout on first task
IF (NPROC > 1) THEN
  IF (MYPROC /= 1) THEN
    OPEN(UNIT=NOUT, FILE='/dev/null')
  ENDIF
ENDIF

! Compute E-W and V-W set sizes
DO JA = INT(SQRT(REAL(NPROC,JPRB))), NPROC
  JB = NPROC / JA
  IF (JA * JB == NPROC) THEN
    NPRGPNS = MAX(JA, JB)
    NPRGPEW = MIN(JA, JB)
    NPRTRW  = MAX(JA, JB)
    NPRTRV  = MIN(JA, JB)
  ENDIF
ENDDO

MYSETV = MOD(MYPROC-1,NPRTRV)+1

! Allocate global spectral arrays
ALLOCATE(ZSPECYG(NFLEVG,NSPEC2G))
ALLOCATE(ZSPECXG(NFLEVG,NSPEC2G))

! Array for storing random perturbations
ALLOCATE(ZRANDSP(NSPEC2G))

! Use a full Gaussian grid
ALLOCATE(NLOEN(NDGL))
NLOEN(:) = 2*NDGL

! Initialise ecTrans
CALL SETUP_TRANS0(KOUT=NOUT, KERR=NERR, KPRINTLEV=0, KMAX_RESOL=1, KPRGPNS=NPRGPNS, &
  &               KPRGPEW=NPRGPEW, KPRTRW=NPRTRW, LDMPOFF=.NOT. LUSE_MPI)
CALL SETUP_TRANS(KSMAX=NSMAX, KDGL=NDGL, KLOEN=NLOEN, LDSPLIT=.TRUE.)

CALL TRANS_INQ(KSPEC2=NSPEC2, KGPTOT=NGPTOT)

! Calculate number of NPROMA blocks
NGPBLKS = (NGPTOT - 1) / NPROMA + 1

! Determine VSET allocation and number of local levels
ALLOCATE(IVSET(NFLEVG))
NFLEV = 0
DO JLEV = 1, NFLEVG
  IVSET(JLEV) = MOD(JLEV,NPRTRV) + 1
  IF (IVSET(JLEV) == MYSETV) THEN
    NFLEV = NFLEV + 1
  ENDIF
ENDDO

! Local spectral arrays
ALLOCATE(ZSPECX(NFLEV,NSPEC2))
ALLOCATE(ZSPECY(NFLEV,NSPEC2))
ALLOCATE(ZSPECP(NFLEV,NSPEC2))

! Temporary grid point array
ALLOCATE(ZGX(NPROMA,NFLEVG,NGPBLKS))

! Prepare perturbations (random numbers between -1 and +1)
IF (MYPROC == 1) THEN
  DO JLEV=1,NFLEVG
    CALL RANDOM_NUMBER(ZRANDSP)
    ZSPECYG(JLEV,:) = (1.0_JPRB-2.0_JPRB*ZRANDSP(:))
    CALL RANDOM_NUMBER(ZRANDSP)
    ZSPECXG(JLEV,:) = (1.0_JPRB-2.0_JPRB*ZRANDSP(:))
  ENDDO
ENDIF

! Distribute global spectral arrays
CALL DIST_SPEC(PSPECG=ZSPECXG, KFDISTG=NFLEVG, KFROM=(/ (1, I = 1, NFLEVG) /), PSPEC=ZSPECX, &
  &            KVSET=IVSET)
CALL DIST_SPEC(PSPECG=ZSPECYG, KFDISTG=NFLEVG, KFROM=(/ (1, I = 1, NFLEVG) /), PSPEC=ZSPECY, &
  &            KVSET=IVSET)

! Calculate DIR_TRANS(INV_TRANS(X))
CALL INV_TRANS(PSPSCALAR=ZSPECX, PGP=ZGX, KPROMA=NPROMA, KVSETSC=IVSET)
CALL DIR_TRANS(PSPSCALAR=ZSPECP, PGP=ZGX, KPROMA=NPROMA, KVSETSC=IVSET)

! Calculate <DIR_TRANS(INV_TRANS(X)), Y>
ZSC1 = SCALPRODSP(ZSPECP, ZSPECY, IVSET, NFLEV, NFLEVG, NSPEC2, NSPEC2G, NSMAX, MYPROC)

ZSPECP = 0.0_JPRB
! Calculate INV_TRANSAD(DIR_TRANSAD(Y))
CALL DIR_TRANSAD(PSPSCALAR=ZSPECY, PGP=ZGX, KPROMA=NPROMA, KVSETSC=IVSET)
CALL INV_TRANSAD(PSPSCALAR=ZSPECP, PGP=ZGX, KPROMA=NPROMA, KVSETSC=IVSET)

! Calculate <X, INV_TRANSAD(DIR_TRANSAD(Y))>
ZSC2 = SCALPRODSP(ZSPECX, ZSPECP, IVSET, NFLEV, NFLEVG, NSPEC2, NSPEC2G, NSMAX, MYPROC)

! If I'm the first task, do the error check
IF (MYPROC == 1) THEN
  ! Calculate relative error between <DIR_TRANS(INV_TRANS(X)), Y> and <X, INV_TRANSAD(DIR_TRANSAD(Y))>
  ZRELATIVE_ERROR = ABS(ZSC1 - ZSC2)/ABS(ZSC1)

  WRITE(NOUT, '(A,1E9.2)') '<Fx,y>  = ', ZSC1
  WRITE(NOUT, '(A,1E9.2)') '<x,F*y> = ', ZSC2
  WRITE(NOUT, '(A,1E9.2)') 'Relative error = ', ZRELATIVE_ERROR

  ! Abort if relative error is > 2000 * machine epsilon
  ! All tested compilers seem to be happy with a threshold of 2000, though it is a bit arbitrary
  IF (ZRELATIVE_ERROR > 2000.0*EPSILON(1.0_JPRB)) THEN
    WRITE(NERR, '(A)') '*******************************'
    WRITE(NERR, '(A)') 'Adjoint test failed'
    WRITE(NERR, '(A)') 'Relative error greater than 2000 * machine epsilon'
    WRITE(NERR, '(1E9.2,A3,1E9.2)') ZRELATIVE_ERROR, ' > ', 2000.0*EPSILON(1.0_JPRB)
    WRITE(NERR, '(A)') '*******************************'
    FLUSH(NERR)
    CALL TRANS_END
    CALL ABORT_TRANS("Adjoint test failed")
  ENDIF
ENDIF

CALL TRANS_END

IF (LUSE_MPI) THEN
  CALL MPL_BARRIER()
  CALL MPL_END
ENDIF

END PROGRAM TEST_ADJOINT
