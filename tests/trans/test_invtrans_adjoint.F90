! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! ==================================================================================================
! INV_TRANS adjoint test
! ==================================================================================================
!
! This program performs a rudimentary check of tangent-linear/adjoint correspondence of 
! INV_TRANS and INV_TRANSAD.
!
! The program checks the correspondence of <INV_TRANS(X1), Y2> and <X1, INV_TRANSAD(Y2)>
! which with infinite precision should match exactly. In practice there is some divergence due to
! rounding errors. In this program we check whether the two expressions are the same to within a
! tolerance of 20000 * machine epsilon.
!
! In this test X1, the "global state vector", is comprised of scalar, vorticity and divergence
! defined on 9 model levels. The correspondence is computed across the whole state vector.
!
! ==================================================================================================

PROGRAM TEST_INVTRANS_ADJOINT

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE, ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC, MPL_BARRIER, MPL_END
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS

IMPLICIT NONE

INTEGER(KIND=JPIM), PARAMETER :: JPTRUNCATION = 159 ! T159 spectral resolution
INTEGER(KIND=JPIM), PARAMETER :: JPPROMA = 16
INTEGER(KIND=JPIM), PARAMETER :: JP_NUMLEVELS_G = 9
INTEGER(KIND=JPIM), PARAMETER :: JPNLAT = 2 * (JPTRUNCATION + 1)

INTEGER(KIND=JPIM) :: INPROC, IMYPROC, IPRGPNS, IPRGPEW, IPRTRW, IPRTRV, IGPTOTG, IGPTOT, IGPBLKS
INTEGER(KIND=JPIM) :: ISPEC2G, ISPEC2
INTEGER(KIND=JPIM) :: INUM_LEVELS
INTEGER(KIND=JPIM) :: IMYSETV, INUMP
INTEGER(KIND=JPIM) :: IVSET(JP_NUMLEVELS_G)
INTEGER(KIND=JPIM), ALLOCATABLE :: MYMS(:)
INTEGER(KIND=JPIM) :: NLOEN(JPNLAT), NASM0(0:JPTRUNCATION)
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

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "inv_trans.h"
#include "inv_transad.h"
#include "gath_grid.h"
#include "dist_grid.h"
#include "gath_spec.h"
#include "dist_spec.h"
#include "trans_end.h"

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
CALL TRANS_INQ(KSPEC2=ISPEC2, KSPEC2G=ISPEC2G, KGPTOT=IGPTOT, KGPTOTG=IGPTOTG, KNUMP=INUMP)
ALLOCATE(MYMS(INUMP))
CALL TRANS_INQ(KMYMS=MYMS, KASM0=NASM0)

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

! ===== Compute invtrans and gather result on proc 1 =====
! i.e. invtrans(rspscalarx, rspvorx, rspdivx) = rgpy

CALL INV_TRANS(PSPSCALAR=ZSPECX, PSPVOR=ZVORX, PSPDIV=ZDIVX, PGP=ZGY, KPROMA=JPPROMA, &
  &            KVSETSC=IVSET, KVSETUV=IVSET)

! ===== Compute: adj_value2 = <invtrans(rspscalarx, rspvorx, rspdivx), rgpx> =====
! i.e. adj_value2 = <rgpy, rgpx>

ADJ_VALUE_1 = SCALPRODGP(ZGY, ZGX)

! ===== Compute adjoint invtrans and gather result on proc 1 =====
! i.e. invtrans_adj(rgpx) = (rspscalary, rspvory, rspdivy)

CALL INV_TRANSAD(PSPSCALAR=ZSPECY, PSPVOR=ZVORY, PSPDIV=ZDIVY, PGP=ZGX, KPROMA=JPPROMA, &
  &              KVSETSC=IVSET, KVSETUV=IVSET)

! ===== Compute: adj_value1 = <(rspscalarx, rspvorx, rspdivx), invtrans_adj(rgpx)> =====
! i.e. adj_value1 = <(rspscalary, rspvory, rspdivy), (rspscalarx, rspvorx, rspdivx)>

ADJ_VALUE_2 = SCALPRODSP(ZSPECX, ZSPECY)
ADJ_VALUE_2 = ADJ_VALUE_2 + SCALPRODSP(ZVORX, ZVORY)
ADJ_VALUE_2 = ADJ_VALUE_2 + SCALPRODSP(ZDIVX, ZDIVY)

! Only task 1 should perform the correctness check
IF (IMYPROC == 1) THEN
  ! ===== Compare inner products =====
  ! i.e. <invtrans_adj(rgpx), (rspscalarx, rspvorx, rspdivx)> == <rgpx, invtrans(rspscalarx, rspvorx, rspdivx)>

  ZRELATIVE_ERROR = ABS(ADJ_VALUE_1 - ADJ_VALUE_2) / ABS(ADJ_VALUE_1)

  WRITE(IOUT, '(A,1E30.15)') '<Fx,y>  = ', ADJ_VALUE_1
  WRITE(IOUT, '(A,1E30.15)') '<x,F*y> = ', ADJ_VALUE_2
  WRITE(IOUT, '(A,1E20.15)') 'Relative error = ', ZRELATIVE_ERROR

  ! Abort if relative error is > 20000 * machine epsilon
  ! All tested compilers seem to be happy with a threshold of 20000, though it is a bit arbitrary
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

CONTAINS

FUNCTION SCALPRODGP(RGP1, RGP2) RESULT(RSC)

  ! Scalar product in spectral space
  REAL(KIND=JPRB) :: RGP1(JPPROMA,3*JP_NUMLEVELS_G,IGPBLKS), RGP2(JPPROMA,3*JP_NUMLEVELS_G,IGPBLKS)
  REAL(KIND=JPRB) :: RSC
  
  INTEGER(KIND=JPIM) :: JLEV, JKGLO, IEND, IBL, JROF
  REAL(KIND=JPRB) :: RGP(JPPROMA,3*JP_NUMLEVELS_G,IGPBLKS), RGPG(IGPTOTG,3*JP_NUMLEVELS_G)
  
  RSC = 0.0_JPRB
  
  !$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JLEV,JKGLO,IEND,IBL,JROF)
  DO JLEV = 1, 3 * JP_NUMLEVELS_G
    DO JKGLO = 1, IGPTOT, JPPROMA
      IEND = MIN(JPPROMA, IGPTOT - JKGLO + 1)
      IBL  = (JKGLO - 1) / JPPROMA+1
      DO JROF = 1, IEND
        RGP(JROF,JLEV,IBL) = RGP1(JROF,JLEV,IBL) * RGP2(JROF,JLEV,IBL)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  
  CALL GATH_GRID(RGPG, JPPROMA, 3*JP_NUMLEVELS_G, ITOGP, PGP=RGP)
  
  IF(IMYPROC == 1) THEN
    RSC = SUM(RGPG)
  ELSE
    RSC = 0.0_JPRB
  ENDIF
  
END FUNCTION SCALPRODGP

FUNCTION SCALPRODSP(PSP1,PSP2) RESULT(PSC)

  ! Scalar product in spectral space
  REAL(KIND=JPRB) :: PSP1(:,:), PSP2(:,:)
  REAL(KIND=JPRB) :: PSC

  INTEGER(KIND=JPIM) :: JMLOC, IM, JIR, JN, INM, JLEV
  REAL(KIND=JPRB) :: ZMFACT,ZSP(INUM_LEVELS,ISPEC2), ZSPG(JP_NUMLEVELS_G,ISPEC2G)

  PSC = 0.0_JPRB
  ZSP(:,:) = 0.0_JPRB

  !$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JLEV,JMLOC,IM,ZMFACT,JIR,JN,INM)
  DO JLEV = 1, INUM_LEVELS
    DO JMLOC = 1, INUMP
      IM = MYMS(JMLOC)
      IF (IM == 0) THEN
        ZMFACT = 1.0_JPRB
      ELSE
        ZMFACT = 2.0_JPRB
      ENDIF
      DO JIR = 0, MIN(1, IM)
        DO JN = IM, JPTRUNCATION
          INM = NASM0(IM) + (JN - IM) * 2 + JIR
          ZSP(JLEV,INM) = PSP1(JLEV,INM) * PSP2(JLEV,INM)*ZMFACT
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  CALL GATH_SPEC(PSPECG=ZSPG, KFGATHG=JP_NUMLEVELS_G, KTO=ITOSP, PSPEC=ZSP, KVSET=IVSET)

  IF (IMYPROC == 1) THEN
    PSC = SUM(ZSPG)
  ELSE
    PSC = 0.0_JPRB
  ENDIF

END FUNCTION SCALPRODSP

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

END PROGRAM TEST_INVTRANS_ADJOINT
