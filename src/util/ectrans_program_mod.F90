module ectrans_program_mod

use iso_fortran_env, only : int32
implicit none
!private
public :: ectrans_program_init, &
          ectrans_program_end, &
          ectrans_print_memory_usage, &
          nthread, nproc, myproc

integer(kind=int32) :: nthread = 1
integer(kind=int32) :: nproc   = 1
integer(kind=int32) :: myproc  = 1

contains

subroutine ectrans_program_init(verbosity, pinning, gstats_config)
  use mpl_module, only: mpl_init, mpl_numproc, mpl_rank, mpl_comm
  use oml_mod ,only : oml_max_threads
  use ec_env_mod, only : ec_putenv
  use yomhook, only: dr_hook_init
  use ectrans_mpi_mod, only: ectrans_mpi_enabled, ectrans_mpi_world_size, ectrans_mpi_world_rank, &
                             ectrans_mpi_world_comm
  use ectrans_device_mod, only : ectrans_device_init, ectrans_device_is_host
  use ectrans_log_mod, only : ectrans_log_init, nout
  use ectrans_gstats_mod, only : ectrans_gstats_init, ectrans_gstats_configuration
  use ectrans_memory_mod, only : allocator

  implicit none
  integer(kind=int32), optional, intent(in) :: verbosity
  logical, optional, intent(in) :: pinning
  type(ectrans_gstats_configuration), intent(in), optional :: gstats_config
  logical :: lpinning
  integer(kind=int32) :: iverbosity
  type(ectrans_gstats_configuration) :: ydgstats_config
#include "abor1.intfb.h"

  if (present(gstats_config)) then
    ydgstats_config = gstats_config
  endif
  call ectrans_device_init()
  iverbosity = 0
  lpinning = .not. ectrans_device_is_host()
  if (present(verbosity)) then
    iverbosity = verbosity
  end if
  if (present(pinning)) then
    lpinning = pinning
  end if

  if (ectrans_mpi_enabled()) then
    call mpl_init(ldinfo=(iverbosity>=1))
    if (ectrans_mpi_world_size() /= mpl_numproc) then
      call abor1fl('ectrans_program_mod.F90',__LINE__,'Mismatch in detected MPI world size between ectrans_mpi_mod and MPL')
    endif
    if (ectrans_mpi_world_rank() /= mpl_rank - 1) then
      call abor1fl('ectrans_program_mod.F90',__LINE__,'Mismatch in detected MPI world rank between ectrans_mpi_mod and MPL')
    endif
  else
    mpl_comm = -1
  endif

  nproc  = ectrans_mpi_world_size()
  myproc = ectrans_mpi_world_rank() + 1

  if (ectrans_mpi_enabled()) then
    block
      use ec_parkind, only : jpim
      use mpl_module, only : mpl_buffer_method, JP_BLOCKING_STANDARD, JP_BLOCKING_BUFFERED
      integer(kind=jpim) :: mp_type  = JP_BLOCKING_BUFFERED ! Message passing type
      integer(kind=jpim) :: mbx_size = 150000000 ! Mailbox size (in bytes) of attached buffer if kmp_type=JP_BLOCKING_BUFFERED
      integer(kind=jpim), allocatable :: nprcids(:)
      integer(kind=jpim) :: jproc
      allocate(nprcids(nproc))
      nprcids = [(jproc, jproc=1,nproc)]
      call mpl_buffer_method(kmp_type=mp_type, kmbx_size=mbx_size, kprocids=nprcids, ldinfo=(iverbosity>=1))
    end block
  endif

  nthread = oml_max_threads()

  call ectrans_log_init(iverbosity)

  !===================================================================================================
  ! Setup allocation strategy
  !===================================================================================================
  if (iverbosity >= 1 .and. ectrans_mpi_world_rank() == 0) then
    call allocator%set_logging(.true.)
    call allocator%set_logging_output_unit(nout)
  endif
  call allocator%set_pinning(lpinning)

  call ectrans_gstats_init(ydgstats_config)

  !call ectrans_default_config%set_default()

  !===================================================================================================
  ! Setup dr_hook
  !===================================================================================================

  call ec_putenv("DR_HOOK_ASSERT_MPI_INITIALIZED=0", overwrite=.true.)
  call dr_hook_init()
end subroutine

subroutine ectrans_program_end()
  use mpl_module, only: mpl_end
  use yomhook, only: dr_hook_end
  use ectrans_mpi_mod, only: ectrans_mpi_enabled
  use ectrans_gstats_mod, only : ectrans_gstats_end
  use ectrans_log_mod, only : ectrans_log_end
  implicit none
  call ectrans_gstats_end()
  call dr_hook_end()
  if (ectrans_mpi_enabled()) then
    call mpl_end(ldmeminfo=.false.)
  endif
  call ectrans_log_end()
end subroutine

subroutine ectrans_print_stack_usage()
  ! NOTE, this seems to be pointless, as within FIAT, the getstackusage function always returns zero.
  !       Consider using the alternative FIAT getstk function.
  use iso_fortran_env, only : int32
  use mpl_module, only: mpl_recv, mpl_send
  use ectrans_mpi_mod, only: ectrans_mpi_world_rank, ectrans_mpi_world_size
  use ectrans_log_mod, only : nout
  implicit none
  integer(kind=int32) :: istack, getstackusage, jproc
  istack = getstackusage()
  if (ectrans_mpi_world_rank() == 0) then
    write(nout,9000) istack
    9000 format("Stack utilisation information",/,&
         &"=============================",/,&
         &"  Task           size(bytes)",/,&
         &"     1",11x,i10)

    do jproc = 2, ectrans_mpi_world_size()
      call mpl_recv(istack, ksource=jproc, ktag=jproc, cdstring='ectrans_print_stack_usage:')
      write(nout,'(i6,11x,i10)') jproc, istack
    enddo
  else
    call mpl_send(istack, kdest=1, ktag=ectrans_mpi_world_rank()+1, cdstring='ectrans_print_stack_usage:')
  endif
end subroutine

subroutine ectrans_print_memory_usage()
  use iso_fortran_env, only : int32
  use mpl_module, only: mpl_barrier, mpl_comm
  use ectrans_mpi_mod, only: ectrans_mpi_enabled
  use ectrans_log_mod, only : nout
  implicit none
#include "ec_meminfo.intfb.h"
  integer :: icomm
  icomm = -1
  if (ectrans_mpi_enabled()) icomm = mpl_comm
  call ec_meminfo(nout, "", kcomm=icomm, kbarr=1, kiotask=-1, kcall=1)
end subroutine

!===================================================================================================

end module