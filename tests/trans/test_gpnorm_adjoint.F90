! (C) Copyright 2024- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! ==================================================================================================
! GPNORM_TRANS adjoint test
! ==================================================================================================
!
! This program performs a rudimentary check of tangent-linear/adjoint correspondence of 
! GPNORM_TRANSTL and GPNORM_TRANSAD.
!
! The program checks the correspondence of <GPNORM_TRANSTL(X1), Y2> and <X1, GPNORM_TRANSAD(Y2)>
! which with infinite precision should match exactly. In practice there is some divergence due to
! rounding errors. In this program we check whether the two expressions are the same to within a
! tolerance of 5000 * machine epsilon.
!
! The check is performed for a grid point array with 10 fields at TCO159 with a block size of 16.
!
! ==================================================================================================

PROGRAM TEST_GPNORM_TRANS_ADJOINT

USE PARKIND1, ONLY: JPIM, JPRB, JPRD
USE MPL_MODULE, ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC, MPL_BARRIER, MPL_END
USE TPM_FIELDS, ONLY: F
USE TPM_GEOMETRY, ONLY: G
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 159 ! T159 spectral resolution
INTEGER(KIND=JPIM), PARAMETER :: JPPROMA = 16
INTEGER(KIND=JPIM), PARAMETER :: JPNUM_FIELDS = 10
INTEGER(KIND=JPIM), PARAMETER :: JPNLAT = 2 * (JPTRUNCATION + 1)

INTEGER(KIND=JPIM) :: INPROC, IMYPROC, IPRGPNS, IPRGPEW, IPRTRW, IPRTRV, IGPTOTG, IGPTOT, IGPBLKS

REAL(KIND=JPRB), ALLOCATABLE :: ZX1(:,:,:), ZX2(:,:,:)
REAL(KIND=JPRB) :: ZY1(JPNUM_FIELDS), ZY2(JPNUM_FIELDS)
REAL(KIND=JPRB), ALLOCATABLE :: ZPRODUCT(:,:,:)
REAL(KIND=JPRB) :: ZPRODUCT_AVE(JPNUM_FIELDS)
INTEGER(KIND=JPIM) :: NLOEN(JPNLAT)

 ! These are not actually used, but they must be passed to GPNORM_TRANS/GPNORM_TRANSAD anyway
REAL(KIND=JPRB) :: ZMIN_DUMMY(JPNUM_FIELDS), ZMAX_DUMMY(JPNUM_FIELDS)

LOGICAL :: LLUSE_MPI
INTEGER(KIND=JPIM) :: IOUT = 6, IERR = 0 ! STDOUT and STDERR
INTEGER(KIND=JPIM) :: JA, JB, JL, JP, JF, JBLK
REAL(KIND=JPRB) :: ZRAND
REAL(KIND=JPRB) :: ZLHS, ZRHS, ZRELATIVE_ERROR
INTEGER(KIND=JPIM) :: N
INTEGER(KIND=JPIM), ALLOCATABLE :: SEED(:)

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "gpnorm_transtl.h"
#include "gpnorm_transad.h"
#include "gpnorm_trans.h"

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

! Compute octahedral latitudes
DO JL = 1, JPNLAT / 2
  NLOEN(JL) = 20 + 4 * (JL - 1)
  NLOEN(JPNLAT - JL + 1) = NLOEN(JL)
END DO

! Initialise ecTrans
CALL SETUP_TRANS0(KOUT=IOUT, KERR=IERR, KPRGPNS=IPRGPNS, KPRGPEW=IPRGPEW, KPRTRW=IPRTRW, &
  &               LDMPOFF=.NOT. LLUSE_MPI)
CALL SETUP_TRANS(KSMAX=JPTRUNCATION, KDGL=2 * (JPTRUNCATION + 1), KLOEN=NLOEN)
CALL TRANS_INQ(KGPTOTG=IGPTOTG, KGPTOT=IGPTOT)

! Initialise grid point arrays
IGPBLKS = (IGPTOT - 1) / JPPROMA + 1
ALLOCATE(ZX1(JPPROMA,JPNUM_FIELDS,IGPBLKS))
ALLOCATE(ZX2(JPPROMA,JPNUM_FIELDS,IGPBLKS))
ALLOCATE(ZPRODUCT(JPPROMA,JPNUM_FIELDS,IGPBLKS))

! Initialise X1 and Y2 with random numbers
DO JP = 1, JPPROMA
  DO JF = 1, JPNUM_FIELDS
    DO JBLK = 1, IGPBLKS
      CALL RANDOM_NUMBER(ZRAND)
      ZX1(JP,JF,JBLK) = (1.0_JPRB - 2.0_JPRB * ZRAND)
    ENDDO
  ENDDO
ENDDO
DO JF = 1, JPNUM_FIELDS
  CALL RANDOM_NUMBER(ZRAND)
  ZY2(JF) = (1.0_JPRB - 2.0_JPRB * ZRAND)
ENDDO

! Calculate TL(X1)
CALL GPNORM_TRANSTL(ZX1, JPNUM_FIELDS, JPPROMA, ZY1)

! Calculate left hand side <TL(X1), Y2>
ZLHS = DOT_PRODUCT(ZY1, ZY2)

! Calculate AD(Y2)
CALL GPNORM_TRANSAD(ZX2, JPNUM_FIELDS, JPPROMA, ZY2)

! Calculate right hand side <X1, AD(Y2)>
ZPRODUCT = ZX1 * ZX2 ! Form the elementwise product in order to calculate the L2 norm of X1 .* X2
F%RW(:) = G%NLOEN(:) ! If we do this, the averaging operation in GPNORM_TRANS becomes an L2 norm
CALL GPNORM_TRANS(ZPRODUCT, JPNUM_FIELDS, JPPROMA, ZPRODUCT_AVE, ZMIN_DUMMY, ZMAX_DUMMY, &
  &               LDAVE_ONLY=.TRUE.)
ZRHS = SUM(ZPRODUCT_AVE) ! Finish the L2 norm over all fields

IF (IMYPROC == 1) THEN
  ! Calculate relative error between LHS and RHS
  ZRELATIVE_ERROR = ABS(ZLHS - ZRHS)/ABS(ZLHS)

  WRITE(IOUT, '(A,1E30.15)') '<TL(X1),Y2>  = ', ZLHS
  WRITE(IOUT, '(A,1E30.15)') '<X1,AD(Y2)> = ', ZRHS
  WRITE(IOUT, '(A,1E20.15)') 'Relative error = ', ZRELATIVE_ERROR

  ! Abort if relative error is > 5000 * machine epsilon
  ! All tested compilers seem to be happy with a threshold of 5000, though it is a bit arbitrary
  IF (ZRELATIVE_ERROR > 5000.0*EPSILON(1.0_JPRB)) THEN
    WRITE(IERR, '(A)') '*******************************'
    WRITE(IERR, '(A)') 'TEST_GPNORM_TRANS_ADJOINT: test failed'
    WRITE(IERR, '(A)') 'Relative error greater than 5000 * machine epsilon'
    WRITE(IERR, '(1E9.2,A3,1E9.2)') ZRELATIVE_ERROR, ' > ', 5000.0*EPSILON(1.0_JPRB)
    WRITE(IERR, '(A)') '*******************************'
    FLUSH(IERR)
    CALL ABORT_TRANS("TEST_GPNORM_TRANS_ADJOINT: test failed")
  ENDIF
ENDIF

IF (LLUSE_MPI) THEN
  CALL MPL_BARRIER()
  CALL MPL_END
ENDIF

CONTAINS

FUNCTION DETECT_MPIRUN() RESULT(LMPI_REQUIRED)
  USE EC_ENV_MOD, ONLY : EC_PUTENV
  LOGICAL :: LMPI_REQUIRED
  INTEGER :: ILEN
  INTEGER, PARAMETER :: NVARS = 4
  CHARACTER(LEN=32), DIMENSION(NVARS) :: CMPIRUN_DETECT
  CHARACTER(LEN=4) :: CLENV
  INTEGER :: IVAR
  
  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  CMPIRUN_DETECT(1) = 'OMPI_COMM_WORLD_SIZE'  ! OPENMPI
  CMPIRUN_DETECT(2) = 'ALPS_APP_PE'           ! CRAY PE
  CMPIRUN_DETECT(3) = 'PMI_SIZE'              ! INTEL
  CMPIRUN_DETECT(4) = 'SLURM_NTASKS'          ! SLURM
  
  LMPI_REQUIRED = .FALSE.
  DO IVAR = 1, NVARS
    CALL GET_ENVIRONMENT_VARIABLE(NAME=TRIM(CMPIRUN_DETECT(IVAR)), LENGTH=ILEN)
    IF (ILEN > 0) THEN
      LMPI_REQUIRED = .TRUE.
      EXIT ! Break
    ENDIF
  ENDDO
  
  CALL GET_ENVIRONMENT_VARIABLE(NAME="ECTRANS_USE_MPI", VALUE=CLENV, LENGTH=ILEN )
  IF (ILEN > 0) THEN
    LMPI_REQUIRED = .TRUE.
    IF( TRIM(CLENV) == "0" .OR. TRIM(CLENV) == "OFF" .OR. TRIM(CLENV) == "OFF" .OR. TRIM(CLENV) == "F" ) THEN
      LMPI_REQUIRED = .FALSE.
    ENDIF
    CALL EC_PUTENV("DR_HOOK_ASSERT_MPI_INITIALIZED=0", OVERWRITE=.TRUE.)
  ENDIF
END FUNCTION

END PROGRAM TEST_GPNORM_TRANS_ADJOINT
