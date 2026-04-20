module ectrans_timer_mod
use iso_fortran_env, only : int32, real64
implicit none
private
public :: ectrans_timer, ectrans_timings

type :: ectrans_timer
  real(kind=real64)   :: start_time
  real(kind=real64)   :: end_time
  integer(kind=int32) :: gstats_id = -1
contains
  procedure :: register_in_gstats => timer_register_in_gstats
  procedure :: start    => timer_start
  procedure :: stop     => timer_stop
  procedure :: elapsed  => timer_elapsed
end type

type :: ectrans_timings
  integer(kind=int32) :: size = 0
  real(kind=real64), allocatable :: times_local(:)
contains
  procedure :: reserve => timings_reserve
  procedure :: push_back => timings_push_back
  procedure :: med => timings_med
  procedure :: min => timings_min
  procedure :: max => timings_max
  procedure :: avg => timings_avg
  procedure :: global_med => timings_global_med
  procedure :: global_min => timings_global_min
  procedure :: global_max => timings_global_max
  procedure :: global_avg => timings_global_avg
end type

contains

function median(vec) result(median_value)
  implicit none
  real(kind=real64), intent(in) :: vec(:)
  real(kind=real64) :: median_value

  real(kind=real64) :: vec_sorted(size(vec))
  real(kind=real64) :: x

  integer(kind=int32) :: i, j, n

  n = size(vec)

  ! Sort in ascending order
  vec_sorted = vec
  do i = 2, n
    x = vec_sorted(i)
    j = i - 1
    do while (j >= 1)
      if (vec_sorted(j) <= x) exit
      vec_sorted(j + 1) = vec_sorted(j)
      j = j - 1
    end do
    vec_sorted(j + 1) = x
  end do

  ! Calculate median according to if there is an even or odd number of elements
  if (mod(n, 2) == 0) then
    median_value = (vec_sorted(n/2) + vec_sorted(n/2+1))/2.0_real64
  else
    median_value = vec_sorted((n+1)/2)
  endif
end function

subroutine timings_reserve(this, ksize)
  class(ectrans_timings), intent(inout) :: this
  integer(kind=int32), intent(in) :: ksize
  real(kind=real64), allocatable :: new_times(:)

  if (allocated(this%times_local)) then
    if (size(this%times_local) < ksize) then
      allocate(new_times(ksize))
      new_times(1:this%size) = this%times_local(1:this%size)
      deallocate(this%times_local)
      allocate(this%times_local(ksize))
      this%times_local = new_times
    end if
  else
    allocate(this%times_local(ksize))
  end if
end subroutine

subroutine timings_push_back(this, time)
  class(ectrans_timings), intent(inout) :: this
  real(kind=real64), intent(in) :: time
  if (.not. allocated(this%times_local)) then
    call this%reserve(1)
  end if
  if (this%size == size(this%times_local)) then
    call this%reserve(max(1, 2*this%size))
  endif
  this%size = this%size + 1
  this%times_local(this%size) = time
end subroutine

function timings_med(this) result(med)
  class(ectrans_timings), intent(in) :: this
  real(kind=real64) :: med
  if (this%size == 0) then
    med = 0.0_real64
  else
    med = median(this%times_local(1:this%size))
  endif
end function

function timings_global_min(this) result(min)
  use mpl_module, only : mpl_allreduce
  use ectrans_mpi_mod, only : ectrans_mpi_enabled
  class(ectrans_timings), intent(inout) :: this
  real(kind=real64) :: min
  if (this%size == 0) then
    min = 0.0_real64
    return
  endif
  min = minval(this%times_local(1:this%size))
  if (ectrans_mpi_enabled()) then
    call mpl_allreduce(min, 'min', ldreprod=.false.)
  endif
end function

function timings_global_max(this) result(max)
  use mpl_module, only : mpl_allreduce
  use ectrans_mpi_mod, only : ectrans_mpi_enabled
  class(ectrans_timings), intent(inout) :: this
  real(kind=real64) :: max
  if (this%size == 0) then
    max = 0.0_real64
    return
  endif
  max = maxval(this%times_local(1:this%size))
  if (ectrans_mpi_enabled()) then
    call mpl_allreduce(max, 'max', ldreprod=.false.)
  endif
end function

function timings_global_avg(this) result(avg)
  use mpl_module, only : mpl_allreduce, mpl_numproc
  use ectrans_mpi_mod, only : ectrans_mpi_enabled
  class(ectrans_timings), intent(inout) :: this
  real(kind=real64) :: avg
  if (this%size == 0) then
    avg = 0.0_real64
    return
  endif
  avg = sum(this%times_local(1:this%size)) / real(this%size, kind=real64)
  if (ectrans_mpi_enabled()) then
    call mpl_allreduce(avg, 'sum', ldreprod=.false.)
    avg = avg / real(mpl_numproc, kind=real64)
  endif
end function

function timings_global_med(this) result(med)
  use mpl_module, only : mpl_allreduce, mpl_numproc
  use ectrans_mpi_mod, only : ectrans_mpi_enabled
  class(ectrans_timings), intent(inout) :: this
  real(kind=real64) :: med
  real(kind=real64), allocatable :: times_global(:)
  if (this%size == 0) then
    med = 0.0_real64
    return
  endif
  if (ectrans_mpi_enabled()) then
    allocate(times_global(this%size))
    times_global = this%times_local(1:this%size)
    call mpl_allreduce(times_global, 'sum', ldreprod=.false.)
    times_global = times_global / real(mpl_numproc, kind=real64)
    med = median(times_global)
    deallocate(times_global)
  else
    med = median(this%times_local(1:this%size))
  endif
end function

function timings_min(this) result(min)
  class(ectrans_timings), intent(in) :: this
  real(kind=real64) :: min
  if (this%size == 0) then
    min = 0.0_real64
  else
    min = minval(this%times_local(1:this%size))
  endif
end function

function timings_max(this) result(max)
  class(ectrans_timings), intent(in) :: this
  real(kind=real64) :: max
  if (this%size == 0) then
    max = 0.0_real64
  else
    max = maxval(this%times_local(1:this%size))
  endif
end function

function timings_avg(this) result(avg)
  class(ectrans_timings), intent(in) :: this
  real(kind=real64) :: avg
  if (this%size == 0) then
    avg = 0.0_real64
  else
    avg = sum(this%times_local(1:this%size)) / real(this%size, kind=real64)
  endif
end function

subroutine timer_register_in_gstats(this,description)
  use ectrans_gstats_mod, only : ectrans_gstats_new_region
  class(ectrans_timer), intent(inout) :: this
  character(len=*), intent(in), optional :: description
  if (present(description)) then
    this%gstats_id = ectrans_gstats_new_region(description=trim(description))
  else
    this%gstats_id = ectrans_gstats_new_region()
  endif
end subroutine

subroutine timer_start(this)
  class(ectrans_timer), intent(inout) :: this
  real(kind=real64), external :: timef ! Timing routine from FIAT returning milliseconds
  if (this%gstats_id >= 0) then
    call gstats(this%gstats_id, 0)
  end if
  this%start_time = timef() / 1000._real64 ! Convert milliseconds to seconds
end subroutine

subroutine timer_stop(this)
  class(ectrans_timer), intent(inout) :: this
  real(kind=real64), external :: timef ! Timing routine from FIAT returning milliseconds
  this%end_time = timef() / 1000._real64 ! Convert milliseconds to seconds
  if (this%gstats_id >= 0) then
    call gstats(this%gstats_id, 1)
  end if
end subroutine

function timer_elapsed(this) result(elapsed)
  class(ectrans_timer), intent(in) :: this
  real(kind=real64) :: elapsed
  elapsed = this%end_time - this%start_time
end function

end module
