module ectrans_log_mod
use iso_fortran_env, only : int32

implicit none
private
public :: nout, nerr, verbosity
public :: ectrans_log_init, ectrans_log_end

integer(kind=int32) :: nerr = 0 ! Unit number for STDERR
integer(kind=int32) :: nout = 6
integer(kind=int32) :: verbosity = 0

contains

subroutine ectrans_log_init(level)
  use ectrans_mpi_mod, only : ectrans_mpi_world_rank
  implicit none
  integer, intent(in), optional :: level
  if (ectrans_mpi_world_rank() > 0) then
    open(unit=nout, file='/dev/null')
  endif
  if (present(level)) then
    verbosity = level
  endif
end subroutine

subroutine ectrans_log_end()
  use ectrans_mpi_mod, only : ectrans_mpi_world_rank
  implicit none
  if (ectrans_mpi_world_rank() > 0) then
    close(unit=nout)
  endif
end subroutine

end module