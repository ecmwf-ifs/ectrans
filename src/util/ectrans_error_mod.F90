module ectrans_error_mod

implicit none
private
public :: ectrans_error_parsing_failed, ectrans_error_print_function
    
abstract interface
  subroutine ectrans_error_print_function(unit)
    integer, intent(in) :: unit
  end subroutine
end interface


interface ectrans_error_parsing_failed
  module procedure ectrans_error_parsing_failed_noprint
  module procedure ectrans_error_parsing_failed_print
end interface


contains

subroutine ectrans_error_abort(errcode)
  ! This routine tries to minimise the error output to avoid overwhelming the user with messages from all MPI ranks.
  ! This is for quick feedback when parsing failed, typically a user error in
  ! the command line arguments or in the configuration file.
  use ectrans_mpi_mod, only : ectrans_mpi_world_rank, ectrans_mpi_world_size
  use mpl_module, only : mpl_abort, mpl_numproc
  implicit none
  integer, intent(in) :: errcode
  external :: ec_sleep
  if (ectrans_mpi_world_size() > 1) then
    if (ectrans_mpi_world_rank() == 0) then
      if (mpl_numproc > 0) then ! if MPL was initialized
        call mpl_abort('ectrans_error_abort')
      endif
      ! Should never reach here since mpl_abort should terminate the program, but just in case.
      error stop 'ectrans_error_abort' 
    else
      ! Give some time for rank 0 to call MPI abort and print the error message before other ranks exit
      ! and potentially cause additional output.
      call ec_sleep(3) ! seconds
      error stop 'ectrans_error_abort'
    end if
  else
    ! In case of only 1 process, just call error stop directly.
    error stop 'ectrans_error_abort'
  endif
end subroutine

subroutine ectrans_error_parsing_failed_noprint(message)
  use ectrans_mpi_mod, only : ectrans_mpi_world_rank
  use ectrans_log_mod, only : nerr
  character(len=*), intent(in) :: message
  if (ectrans_mpi_world_rank() == 0) then
    write(nerr,"(a)") trim(message)
  endif
  call ectrans_error_abort(1)
end subroutine

subroutine ectrans_error_parsing_failed_print(message, print_function)
  use ectrans_mpi_mod, only : ectrans_mpi_world_rank
  use ectrans_log_mod, only : nerr
  character(len=*), intent(in) :: message
  procedure(ectrans_error_print_function) :: print_function
  if (ectrans_mpi_world_rank() == 0) then
    write(nerr,"(a)") trim(message)
    write(nerr,"(a)") "--------------------------------------------------------------------------"
    call print_function(unit=nerr)
    write(nerr,"(a)") "--------------------------------------------------------------------------"
  endif
  call ectrans_error_abort(1)
end subroutine

end module