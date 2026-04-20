module ectrans_checksum_mod
  use iso_c_binding, only : c_ptr, c_null_ptr

implicit none
private

type, public :: ectrans_checksum_file_writer
  character(len=:), allocatable :: file_name
  integer :: file_unit
  logical :: file_open = .false.
contains
  procedure, public :: open  => open_file
  procedure, public :: close => close_file
  procedure, public :: write_iteration_separator
  procedure, public :: write_checksums_pgp
  procedure, public :: write_checksums_pgp_uv_3a_2
  procedure, public :: write_checksums_psp
  procedure, public :: write_checksums_psp_3a_2
  final :: ectrans_checksum_file_writer_final
end type ectrans_checksum_file_writer

type, private :: memory_buffer
  type(c_ptr) :: data = c_null_ptr
  integer :: size = 0
contains
  procedure :: reserve => memory_buffer_reserve
end type memory_buffer

! A growable memory buffer for reuse between checksum calculations
! to avoid repeated malloc/free calls which can be expensive.
type(memory_buffer), private :: buffer
interface get_buffer
  module procedure get_buffer_real32
  module procedure get_buffer_real64
end interface


contains

!===================================================================================================

subroutine deallocate_buffer_ptr(ptr)
  use iso_c_binding, only : c_ptr, c_associated, c_null_ptr
  implicit none
  type(c_ptr), intent(inout) :: ptr
  interface
    subroutine c_free(ptr) bind(C, name="free")
      import :: c_ptr
      type(c_ptr), value :: ptr
    end subroutine c_free
  end interface
  if (c_associated(ptr)) then
    call c_free(ptr)
    ptr = c_null_ptr
  endif
end subroutine deallocate_buffer_ptr

subroutine allocate_buffer_ptr(ptr, bytes)
  use iso_c_binding, only : c_ptr, c_size_t
  implicit none
  type(c_ptr), intent(inout) :: ptr
  integer(c_size_t), intent(in) :: bytes
  interface
    function c_malloc(size) bind(C, name="malloc")
      import :: c_ptr, c_size_t
      type(c_ptr) :: c_malloc
      integer(c_size_t), value :: size
    end function c_malloc
  end interface
  ptr = c_malloc(bytes)
end subroutine allocate_buffer_ptr

subroutine memory_buffer_reserve(this, bytes)
  use iso_c_binding, only : c_size_t
  implicit none
  class(memory_buffer), intent(inout) :: this
  integer(c_size_t), intent(in) :: bytes
  if (bytes > this%size) then
    call deallocate_buffer_ptr(this%data)
    call allocate_buffer_ptr(this%data, bytes)
    this%size = bytes
  endif
end subroutine memory_buffer_reserve

subroutine get_buffer_real32(array, shape)
  use iso_c_binding, only : c_size_t, c_double, c_f_pointer
  use iso_fortran_env, only : real32
  implicit none
  real(real32), pointer, intent(inout) :: array(:,:)
  integer, intent(in) :: shape(2)
  integer(c_size_t) :: bytes
  bytes = storage_size(array, c_size_t) / 8_c_size_t * int(shape(1),c_size_t) * int(shape(2),c_size_t)
  call buffer%reserve(bytes)
  call c_f_pointer(buffer%data, array, shape)
end subroutine get_buffer_real32

subroutine get_buffer_real64(array, shape)
  use iso_c_binding, only : c_size_t, c_double, c_f_pointer
  use iso_fortran_env, only : real64
  implicit none
  real(real64), pointer, intent(inout) :: array(:,:)
  integer, intent(in) :: shape(2)
  integer(c_size_t) :: bytes
  bytes = storage_size(array, c_size_t) / 8_c_size_t * int(shape(1),c_size_t) * int(shape(2),c_size_t)
  call buffer%reserve(bytes)
  call c_f_pointer(buffer%data, array, shape)
end subroutine get_buffer_real64

subroutine assert_file_open(this)
#include "abor1.intfb.h"
  class(ectrans_checksum_file_writer), intent(in) :: this
  character(len=1024) :: errmsg
  if (.not. this%file_open) then
    write(errmsg,'(A,A,A,I0)') "File not open: ", this%file_name, ' unit=',this%file_unit
    call abor1(errmsg)
  endif
end subroutine assert_file_open

subroutine trans_inq_myproc(myproc)
  use ectrans_mpi_mod, only : ectrans_mpi_enabled
  use mpl_module, only : mpl_rank
  integer, intent(out) :: myproc
  ! In the future we should be able to trans_inq this.
  ! For now ectrans uses MPL communicators

  myproc = 1
  if (ectrans_mpi_enabled()) then
    myproc = mpl_rank
  endif
end subroutine trans_inq_myproc

subroutine open_file(this, filename, append)
  use ectrans_log_mod, only : nout, verbosity
  class(ectrans_checksum_file_writer), intent(inout) :: this
  character(len=*), intent(in) :: filename
  logical, optional :: append
  logical :: lappend
  integer :: myproc
  character(len=1024) :: errmsg
  call trans_inq_myproc(myproc)
  if (myproc /= 1) return
  if (this%file_open) then
    write(errmsg,'(A,A,A,I0)') "File already open: ", this%file_name, ' unit=',this%file_unit
    call abor1(errmsg)
  endif
  this%file_name = trim(filename)
  lappend = .false.
  if (present(append)) then
    lappend = append
  endif
  if (lappend) inquire(file=this%file_name, exist=lappend)
  if (lappend) then
    open(newunit=this%file_unit, file=this%file_name, status="old", position="append", action="write")
    if (verbosity >= 1) then
      write(nout,'(A,A,A,I0)') "ectrans_checksum_file_writer opened for append : ", this%file_name, ' unit=', this%file_unit
    endif
    write(nout,'(A,A,A,I0)') "ectrans_checksum_file_writer re-opened for append : ", this%file_name, ' unit=', this%file_unit
  else
    open(newunit=this%file_unit, file=this%file_name, action="write")
    if (verbosity >= 1) then
      write(nout,'(A,A,A,I0)') "ectrans_checksum_file_writer opened for write : ", this%file_name, ' unit=', this%file_unit
    endif
  endif
  this%file_open = .true.
end subroutine open_file

subroutine close_file(this)
  use ectrans_log_mod, only : nout, verbosity
  class(ectrans_checksum_file_writer), intent(inout) :: this
  integer :: myproc
  call trans_inq_myproc(myproc)
  if (myproc /= 1) return
  if (this%file_open) then
    close(this%file_unit)
    if (verbosity >= 1) then
      write(nout,'(A,A,A,I0)') "ectrans_checksum_file_writer closed file ", this%file_name, ' unit=', this%file_unit
    endif
  endif
  this%file_open = .false.
end subroutine close_file

subroutine ectrans_checksum_file_writer_final(this)
  type(ectrans_checksum_file_writer), intent(inout) :: this
  call this%close()
end subroutine ectrans_checksum_file_writer_final


subroutine write_iteration_separator(this,iteration)
  use ec_parkind, only : jpim
  class(ectrans_checksum_file_writer), intent(inout) :: this
  integer, intent(in) :: iteration
  integer(kind=jpim) :: myproc
  call trans_inq_myproc(myproc)
  if (myproc == 1) then
    call assert_file_open(this)
    write(this%file_unit,'(a)') "--------------------------------------------"
    write(this%file_unit,'(a, i0)') "Iteration ", iteration
    write(this%file_unit,'(a)') "--------------------------------------------"
  endif
end subroutine write_iteration_separator

subroutine write_checksums_pgp(this, zgp)
  use parkind1, only : jprb, jpib, jpim
  use ectrans_log_mod, only : nout
  implicit none
  class(ectrans_checksum_file_writer), intent(inout) :: this
  real(kind=jprb), intent(in) :: zgp(:,:,:)
  integer(kind=jpib) :: icrc
  integer(kind=jpim) :: jfld
  real(kind=jprb), pointer :: gfld(:,:)

  integer(kind=jpim) :: nproma   ! size of nproma
  integer(kind=jpim) :: ngptotg
  integer(kind=jpim) :: myproc

#include "gath_grid.h"
#include "trans_inq.h"

  call trans_inq(kgptotg=ngptotg)
  nproma = size(zgp, 1)
  call trans_inq_myproc(myproc)

  if (myproc == 1) then
    call get_buffer(gfld, [ngptotg, 1])
  endif

  icrc = 0
  do jfld = 1, size(zgp, 2)
    call gath_grid(pgpg=gfld, kproma=nproma, kfgathg=1, kto=(/1/), kresol=1, &
      &            pgp=zgp(:,jfld:jfld,:))
    if (myproc == 1) then
      call crc64(gfld(:,:), int(size(gfld(:,:)) * kind(gfld), 8), icrc)
      write(this%file_unit, '(a," (",i0,") = ",z16.16)') "zgp", jfld, icrc
    endif
  enddo

end subroutine write_checksums_pgp

!===================================================================================================

subroutine write_checksums_pgp_uv_3a_2(this, zgpuv, zgp3a, zgp2)
  use parkind1, only : jprb, jpib, jpim
  use ectrans_log_mod, only : nout
  implicit none
  class(ectrans_checksum_file_writer), intent(inout) :: this
  real(kind=jprb), intent(in) :: zgpuv(:,:,:,:)
  real(kind=jprb), intent(in) :: zgp3a(:,:,:,:)
  real(kind=jprb), intent(in) :: zgp2(:,:,:)

  integer(kind=jpib) :: icrc
  integer(kind=jpim) :: jlev, jfld
  real(kind=jprb), pointer :: gfld(:,:)

  integer(kind=jpim) :: myproc
  integer(kind=jpim) :: nproma   ! size of nproma
  integer(kind=jpim) :: ngptotg

#include "gath_grid.h"
#include "trans_inq.h"

  call trans_inq_myproc(myproc)
  nproma = size(zgpuv, 1)
  call trans_inq(kgptotg=ngptotg)

  if (myproc == 1) then
    call assert_file_open(this)
    call get_buffer(gfld, [ngptotg, 1])
  endif

  icrc = 0
  do jfld = 1, size(zgpuv, 3)
    do jlev = 1, size(zgpuv, 2)
      call gath_grid(pgpg=gfld, kproma=nproma, kfgathg=1, kto=(/1/), kresol=1, &
        &            pgp=zgpuv(:,jlev:jlev,jfld,:))
      if (myproc == 1) then
        call crc64(gfld(:,:), int(size(gfld(:,:)) * kind(gfld), 8), icrc)
        write(this%file_unit, '(a," (",i0,", ",i0,") = ",z16.16)') "zgpuv", jlev, jfld, icrc
      endif
    enddo
  enddo

  icrc = 0
  do jfld = 1, size(zgp3a, 3)
    do jlev = 1, size(zgp3a, 2)
      call gath_grid(pgpg=gfld, kproma=nproma, kfgathg=1, kto=(/1/), kresol=1, &
        &            pgp=zgp3a(:,jlev:jlev,jfld,:))
      if (myproc == 1) then
        call crc64(gfld(:,:), int(size(gfld(:,:)) * kind(gfld), 8), icrc)
        write(this%file_unit, '(a," (",i0,", ",i0,") = ",z16.16)') "zgp3a", jlev, jfld, icrc
      endif
    enddo
  enddo

  icrc = 0
  do jfld = 1, size(zgp2, 2)
    call gath_grid(pgpg=gfld, kproma=nproma, kfgathg=1, kto=(/1/), kresol=1, &
      &            pgp=zgp2(:,jfld:jfld,:))
    if (myproc == 1) then
      call crc64(gfld(:,:), int(size(gfld(:,:)) * kind(gfld), 8), icrc)
      write(this%file_unit, '(a," (",i0,") = ",z16.16)') "zgp2", jfld, icrc
    endif
  enddo

end subroutine write_checksums_pgp_uv_3a_2

!===================================================================================================

subroutine write_checksums_psp(this, ivset, ivsetsc, zspvor, zspdiv, zspscalar)
  use parkind1, only : jprb, jpib, jpim
  use ectrans_log_mod, only : nout
  implicit none
  class(ectrans_checksum_file_writer), intent(inout) :: this

  integer(kind=jpim), intent(in) :: ivset(:)
  integer(kind=jpim), intent(in) :: ivsetsc(:)
  real(kind=jprb), intent(in) :: zspvor(:,:)
  real(kind=jprb), intent(in) :: zspdiv(:,:)
  real(kind=jprb), intent(in) :: zspscalar(:,:)
  integer(kind=jpim) :: numfld, i
  integer(kind=jpib) :: icrc
  real(kind=jprb), pointer :: gspfld(:,:)

  integer(kind=jpim) :: myproc
  integer(kind=jpim) :: nspec2g

#include "gath_spec.h"
#include "trans_inq.h"

  call trans_inq_myproc(myproc)
  call trans_inq(kspec2g=nspec2g)

  if (myproc == 1) then
    call assert_file_open(this)
    call get_buffer(gspfld, [max(size(ivset), size(ivsetsc)), nspec2g])
  endif

  numfld = size(ivset)
  if (myproc == 1) then
    call gath_spec(pspecg=gspfld(1:numfld,:), kfgathg=numfld, kto=[(1, i = 1, numfld)], &
      &            kvset=ivset, pspec=zspvor)
    icrc = 0
    call crc64(gspfld(1:numfld,:), int(size(gspfld(1:numfld,:)) * kind(gspfld), 8), icrc)
    write(this%file_unit, '(a," = ",z16.16)') "zspvor", icrc
  else
    call gath_spec(kfgathg=numfld, kto=[(1, i = 1, numfld)], kvset=ivset, pspec=zspvor)
  endif

  if (myproc == 1) then
    call gath_spec(pspecg=gspfld(1:numfld,:), kfgathg=numfld, kto=[(1, i = 1, numfld)], &
      &            kvset=ivset, pspec=zspdiv)
    icrc = 0
    call crc64(gspfld(1:numfld,:), int(size(gspfld(1:numfld,:)) * kind(gspfld), 8), icrc)
    write(this%file_unit, '(a," = ",z16.16)') "zspdiv", icrc
  else
    call gath_spec(kfgathg=numfld, kto=[(1, i = 1, numfld)], kvset=ivset, pspec=zspdiv)
  endif

  numfld = size(ivsetsc)
  if (myproc == 1) then
    call gath_spec(pspecg=gspfld(1:numfld,:), kfgathg=numfld, kto=[(1, i = 1, numfld)], &
      &            kvset=ivsetsc, pspec=zspscalar)
    icrc = 0
    call crc64(gspfld(1:numfld,:), int(size(gspfld(1:numfld,:)) * kind(gspfld), 8), icrc)
    write(this%file_unit, '(a," = ",z16.16)') "zspscalar", icrc
  else
    call gath_spec(kfgathg=numfld, kto=[(1, i = 1, numfld)], kvset=ivsetsc, pspec=zspscalar)
  endif

end subroutine write_checksums_psp

!===================================================================================================

subroutine write_checksums_psp_3a_2(this,  &
                        & ivset, ivsetsc2,              &
                        & zspvor, zspdiv,               &
                        & zspsc3a, zspsc2)
  use parkind1, only : jprb, jpib, jpim
  use ectrans_log_mod, only : nout
  implicit none
  class(ectrans_checksum_file_writer), intent(inout) :: this
  integer(kind=jpim), intent(in) :: ivset(:)
  integer(kind=jpim), intent(in) :: ivsetsc2(:)
  real(kind=jprb), intent(in) :: zspvor(:,:)
  real(kind=jprb), intent(in) :: zspdiv(:,:)
  real(kind=jprb), intent(in) :: zspsc3a(:,:,:)
  real(kind=jprb), intent(in) :: zspsc2(:,:)

  integer(kind=jpim) :: numfld, jfld, i
  integer(kind=jpib) :: icrc
  integer(kind=jpim) :: myproc
  integer(kind=jpim) :: nspec2g

  real(kind=jprb), pointer :: gspfld(:,:)
#include "gath_spec.h"
#include "trans_inq.h"

  call trans_inq_myproc(myproc)
  call trans_inq(kspec2g=nspec2g)

  if (myproc == 1) then
    call assert_file_open(this)
    call get_buffer(gspfld, [max(size(ivset), 1), nspec2g]) ! size(ivsetsc2) is always 1
  endif

  numfld = size(ivset)
  if (myproc == 1) then
    call gath_spec(pspecg=gspfld(1:numfld,:), kfgathg=numfld, kto=[(1, i = 1, numfld)], &
      &            kvset=ivset, pspec=zspvor)
    icrc = 0
    call crc64(gspfld(1:numfld,:), int(size(gspfld(1:numfld,:)) * kind(gspfld), 8), icrc)
    write(this%file_unit, '(a," = ",z16.16)') "zspvor", icrc
  else
    call gath_spec(kfgathg=numfld, kto=[(1, i = 1, numfld)], kvset=ivset, pspec=zspvor)
  endif

  if (myproc == 1) then
    call gath_spec(pspecg=gspfld(1:numfld,:), kfgathg=numfld, kto=[(1, i = 1, numfld)], &
      &            kvset=ivset, pspec=zspdiv)
    icrc = 0
    call crc64(gspfld(1:numfld,:), int(size(gspfld(1:numfld,:)) * kind(gspfld), 8), icrc)
    write(this%file_unit, '(a," = ",z16.16)') "zspdiv", icrc
  else
    call gath_spec(kfgathg=numfld, kto=[(1, i = 1, numfld)], kvset=ivset, pspec=zspdiv)
  endif

  do jfld = 1, size(zspsc3a, 3)
    if (myproc == 1) then
      call gath_spec(pspecg=gspfld(1:numfld,:), kfgathg=numfld, kto=[(1, i = 1, numfld)], &
        &            kvset=ivset, pspec=zspsc3a(:,:,jfld))
      icrc = 0
      call crc64(gspfld(1:numfld,:), int(size(gspfld(1:numfld,:)) * kind(gspfld), 8), icrc)
      write(this%file_unit, '(a,"(",i0,") = ",z16.16)') "zspsc3a", jfld, icrc
    else
      call gath_spec(kfgathg=numfld, kto=[(1, i = 1, numfld)], kvset=ivset, pspec=zspsc3a(:,:,jfld))
    endif
  enddo

  if (myproc == 1) then
    call gath_spec(pspecg=gspfld(1:1,:), kfgathg=1, kto=[1], kvset=ivsetsc2, pspec=zspsc2)
    icrc = 0
    call crc64(gspfld(1,:), int(size(gspfld(1,:)) * kind(gspfld), 8), icrc)
    write(this%file_unit, '(a," = ",z16.16)') "zspsc2", icrc
  else
    call gath_spec(kfgathg=1, kto=[1], kvset=ivsetsc2, pspec=zspsc2)
  endif

end subroutine write_checksums_psp_3a_2

end module