! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM TEST_SETUP_TRANS0
USE EC_PARKIND  ,ONLY : JPIM
USE MPL_MODULE  ,ONLY : MPL_INIT, MPL_END, MPL_BARRIER, MPL_MYRANK, MPL_NPROC
USE ABORT_TRANS_MOD, ONLY : ABORT_TRANS
USE YOMHOOK, ONLY : JPHOOK, DR_HOOK

IMPLICIT NONE

INTEGER(KIND=JPIM) :: NPROC,NPRGPNS,NPRGPEW,NPRTRW,NPRTRV
INTEGER(KIND=JPIM) :: NOUT, MYPROC
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
CHARACTER(LEN=6) :: CLNAME
LOGICAL :: LUSE_MPI

#include "setup_trans0.h"

LUSE_MPI = detect_mpirun()

IF(LUSE_MPI) THEN
  CALL MPL_INIT
  MYPROC = MPL_MYRANK()
  NPROC = MPL_NPROC()
  NOUT = 20
  WRITE(CLNAME,'(A,I2.2)') 'OUT.',MYPROC
  OPEN(NOUT,FILE=CLNAME)
ELSE
  NOUT   = 6
  MYPROC = 1
  NPROC  = 1
ENDIF

CALL DR_HOOK('PROGRAM', 0, ZHOOK_HANDLE)

! ======================================================================================
! NPROC must match NPRGPNS * NPRGPEW
NPRTRV  = 1
NPRTRW  = NPROC / NPRTRV
NPRGPEW = 1
NPRGPNS = NPROC
! ======================================================================================

IF (MYPROC == 1) WRITE(NOUT,*) ' LUSE_MPI= ',LUSE_MPI

IF(NPROC /= NPRTRW*NPRTRV) THEN
  PRINT *,'NPRGPNS,NPRGPEW,NPRTRW,NPRTRV ',NPRGPNS,NPRGPEW,NPRTRW,NPRTRV
  CALL ABORT_TRANS('NPRGPNS*NPRGPEW /= NPRTRW*NPRTRV')
ENDIF

CALL SETUP_TRANS0(KOUT=NOUT,KERR=0,KPRINTLEV=2, &
 & KMAX_RESOL=1,&
 & LDEQ_REGIONS=.TRUE., &
 & KPRGPNS=NPRGPNS, KPRGPEW=NPRGPEW, KPRTRW=NPRTRW,&
 & LDMPOFF=.NOT.LUSE_MPI)

CALL DR_HOOK('PROGRAM', 1, ZHOOK_HANDLE)

IF(LUSE_MPI) THEN
 CALL MPL_BARRIER()
 CALL MPL_END
ENDIF

CONTAINS

function detect_mpirun() result(lmpi_required)
  use ec_env_mod, only : ec_putenv
  logical :: lmpi_required
  integer :: ilen
  integer, parameter :: nvars = 4
  character(len=32), dimension(nvars) :: cmpirun_detect
  character(len=4) :: clenv
  integer :: ivar

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  cmpirun_detect(1) = 'OMPI_COMM_WORLD_SIZE'  ! openmpi
  cmpirun_detect(2) = 'ALPS_APP_PE'           ! cray pe
  cmpirun_detect(3) = 'PMI_SIZE'              ! intel
  cmpirun_detect(4) = 'SLURM_NTASKS'          ! slurm

  lmpi_required = .false.
  do ivar = 1, nvars
    call get_environment_variable(name=trim(cmpirun_detect(ivar)), length=ilen)
    if (ilen > 0) then
      lmpi_required = .true.
      exit ! break
    endif
  enddo

  call get_environment_variable(name="ECTRANS_USE_MPI", value=clenv, length=ilen )
  if (ilen > 0) then
      lmpi_required = .true.
      if( trim(clenv) == "0" .or. trim(clenv) == "OFF" .or. trim(CLENV) == "off" .or. trim(clenv) == "F" ) then
        lmpi_required = .false.
      endif
      call ec_putenv("DR_HOOK_ASSERT_MPI_INITIALIZED=0", overwrite=.true.)
  endif
end function


END PROGRAM
