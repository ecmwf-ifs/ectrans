module ectrans_command_line_parser_mod

implicit none
private
public :: ectrans_command_line_parser
    
abstract interface
  subroutine print_help_signature(unit)
    integer, intent(in) :: unit
  end subroutine
end interface

type :: ectrans_command_line_parser
  integer :: iarg = 0
  integer :: narg = 0
  character(len=128) :: carg
  procedure(print_help_signature), nopass, pointer :: print_help_ptr => null()

contains
  procedure :: setup          => command_line_parser_setup
  procedure :: next_arg       => command_line_parser_next_arg
  procedure :: get_int_value  => command_line_parser_get_value_int
  procedure :: get_str_value  => command_line_parser_get_value_str
  procedure :: parsing_failed => command_line_parser_parsing_failed
  procedure :: print_help     => command_line_parser_print_help
end type

contains

subroutine command_line_parser_parsing_failed(this, message)
  use ectrans_error_mod, only : ectrans_error_parsing_failed
  class(ectrans_command_line_parser), intent(inout) :: this
  character(len=*), intent(in) :: message
  call ectrans_error_parsing_failed(message, this%print_help_ptr)
end subroutine

subroutine command_line_parser_setup(this, print_help)
  implicit none
  class(ectrans_command_line_parser), intent(inout) :: this
  procedure(print_help_signature) :: print_help
  this%print_help_ptr => print_help
  this%narg = command_argument_count()
end subroutine command_line_parser_setup

function command_line_parser_next_arg(this, next_arg) result(has_next)
  implicit none
  class(ectrans_command_line_parser), intent(inout) :: this
  character(len=128) :: next_arg
  logical :: has_next
  if (this%iarg < this%narg) then
    has_next = .true.
    this%iarg = this%iarg + 1
    call get_command_argument(this%iarg, this%carg)
    next_arg = this%carg
  else
    next_arg = ""
    has_next = .false.
  end if
end function

subroutine command_line_parser_print_help(this, unit)
  use ectrans_mpi_mod, only : ectrans_mpi_world_size, ectrans_mpi_world_rank
  use ectrans_log_mod, only : nout
  use mpl_module, only : mpl_end, mpl_numproc
  implicit none
  class(ectrans_command_line_parser), intent(in) :: this
  integer, optional, intent(in) :: unit
  if (ectrans_mpi_world_rank() == 0) then
    if (present(unit)) then
      call this%print_help_ptr(unit)
    else
      call this%print_help_ptr(nout)
    end if
  end if
  if (ectrans_mpi_world_size() > 1) then
    if (mpl_numproc > 0) then
      call mpl_end()
    end if
  endif
end subroutine command_line_parser_print_help

subroutine str2int(str, int, stat)
  implicit none
  character(len=*), intent(in) :: str
  integer, intent(out) :: int
  integer, intent(out) :: stat
  read(str, *, iostat=stat) int
end subroutine str2int

function command_line_parser_get_value_str(this) result(value)
  implicit none
  class(ectrans_command_line_parser), intent(inout) :: this
  character(len=128) :: name
  character(len=128) :: value
  name = this%carg
  this%iarg = this%iarg + 1
  call get_command_argument(this%iarg, value)
  if (value == "") then
    call this%parsing_failed("Invalid argument for " // trim(name) // ": no value provided")
  end if
end function

function command_line_parser_get_value_int(this) result(value)
  implicit none
  class(ectrans_command_line_parser), intent(inout) :: this
  integer :: value
  character(len=128) :: name
  character(len=128) :: value_str

  integer :: stat
  name = this%carg

  this%iarg = this%iarg + 1
  call get_command_argument(this%iarg, value_str)
  if (value_str == "") then
    call this%parsing_failed("Invalid argument for " // trim(name) // ": no value provided")
  end if

  call str2int(value_str, value, stat)

  if (stat /= 0) then
    call this%parsing_failed("Invalid argument for " // trim(name) // ": " // trim(value_str))
  end if
end function

end module