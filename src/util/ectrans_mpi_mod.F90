module ectrans_mpi_mod
use iso_fortran_env, only : int32

implicit none

private

public :: ectrans_mpi_detect, ectrans_mpi_enabled, ectrans_mpi_world_size, ectrans_mpi_world_rank, &
          ectrans_mpi_world_comm

logical,             save :: ectrans_mpi_enabled_   = .false.
integer(kind=int32), save :: ectrans_mpi_world_size_ = 1
integer(kind=int32), save :: ectrans_mpi_world_rank_ = 0
integer(kind=int32), save :: ectrans_mpi_world_comm_ = -1

contains

subroutine ectrans_mpi_detect()
  ! Detect if MPI is enabled, and set values for world size, rank and communicator, without MPI being required to be initialized.
  ! The detection of MPI is based on the presence of environment variables that
  ! are typically set when using mpirun, srun, aprun, ... to launch the program.
  ! Additionally, the ECTRANS_USE_MPI environment variable can be used to force enable or disable MPI support.

  use mpl_mpif, only : mpi_comm_world
  implicit none
  integer(kind=int32) :: ilen
  integer(kind=int32), parameter :: nvars = 4
  character(len=32), dimension(nvars) :: cmpirun_detect
  character(len=4) :: clenv
  integer(kind=int32) :: ivar
  integer(kind=int32), external :: ec_mpirank
  integer(kind=int32), external :: ec_mpisize
  logical, save :: already_called = .false.
  if (already_called) return
  already_called = .true.

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  cmpirun_detect(1) = 'OMPI_COMM_WORLD_SIZE'  ! openmpi
  cmpirun_detect(2) = 'ALPS_APP_PE'           ! cray pe
  cmpirun_detect(3) = 'PMI_SIZE'              ! intel
  cmpirun_detect(4) = 'SLURM_STEP_NUM_TASKS'  ! slurm

  ectrans_mpi_enabled_ = .false.
  do ivar = 1, nvars
    call get_environment_variable(name=trim(cmpirun_detect(ivar)), length=ilen)
    if (ilen > 0) then
      ectrans_mpi_enabled_ = .true.
      exit ! break
    endif
  enddo

  call get_environment_variable(name="ECTRANS_USE_MPI", value=clenv, length=ilen )
  if (ilen > 0) then
    ectrans_mpi_enabled_ = .true.
    if( trim(clenv) == "0" .or. trim(clenv) == "OFF" .or. trim(clenv) == "off" .or. trim(clenv) == "F" ) then
      ectrans_mpi_enabled_ = .false.
    endif
  endif

  if (ectrans_mpi_enabled_) then
    ectrans_mpi_world_size_ = ec_mpisize()
    ectrans_mpi_world_rank_ = ec_mpirank()
    ectrans_mpi_world_comm_ = mpi_comm_world
  end if
end subroutine

function ectrans_mpi_enabled() result(enabled)
  logical :: enabled
  call ectrans_mpi_detect()
  enabled = ectrans_mpi_enabled_
end function

function ectrans_mpi_world_size() result(world_size)
  integer(kind=int32) :: world_size
  call ectrans_mpi_detect()
  world_size = ectrans_mpi_world_size_
end function

function ectrans_mpi_world_rank() result(world_rank)
  integer(kind=int32) :: world_rank
  call ectrans_mpi_detect()
  world_rank = ectrans_mpi_world_rank_
end function

function ectrans_mpi_world_comm() result(world_comm)
  integer(kind=int32) :: world_comm
  call ectrans_mpi_detect()
  world_comm = ectrans_mpi_world_comm_
end function

end module