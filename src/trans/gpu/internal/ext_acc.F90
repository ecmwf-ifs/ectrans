! (C) Copyright 2022- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
module openacc_ext_type
  use iso_c_binding
  implicit none
  private
  public :: ext_acc_arr_desc

  ! to my knowledge, this cannot be part of openacc_ext
  type ext_acc_arr_desc
    integer(c_size_t) :: ptr, sz
  end type
end module
module openacc_ext
  use iso_c_binding
  use openacc
  use openacc_ext_type
  implicit none

  private
  public :: ext_acc_pass, ext_acc_create, ext_acc_copyin, ext_acc_copyout, &
    & ext_acc_delete, ext_acc_arr_desc, acc_handle_kind

  type common_pointer_descr
    type(c_ptr) :: ptr
    integer(c_size_t) :: sz
  end type

  interface ext_acc_pass
    function ext_acc_pass_2d_r4(arr) result(ret)
      use openacc_ext_type
      implicit none
      type(ext_acc_arr_desc) :: ret
      real(4), intent(in) :: arr(:,:)
    end function
    function ext_acc_pass_3d_r4(arr) result(ret)
      use openacc_ext_type
      implicit none
      type(ext_acc_arr_desc) :: ret
      real(4), intent(in) :: arr(:,:,:)
    end function
    function ext_acc_pass_4d_r4(arr) result(ret)
      use openacc_ext_type
      implicit none
      type(ext_acc_arr_desc) :: ret
      real(4), intent(in) :: arr(:,:,:,:)
    end function
    function ext_acc_pass_2d_r8(arr) result(ret)
      use openacc_ext_type
      implicit none
      type(ext_acc_arr_desc) :: ret
      real(8), intent(in) :: arr(:,:)
    end function
    function ext_acc_pass_3d_r8(arr) result(ret)
      use openacc_ext_type
      implicit none
      type(ext_acc_arr_desc) :: ret
      real(8), intent(in) :: arr(:,:,:)
    end function
    function ext_acc_pass_4d_r8(arr) result(ret)
      use openacc_ext_type
      implicit none
      type(ext_acc_arr_desc) :: ret
      real(8), intent(in) :: arr(:,:,:,:)
    end function
  end interface
contains

  function ext_acc_pass_2d_r4(arr) result(ret)
    implicit none
    type(ext_acc_arr_desc) :: ret
    real(4), intent(in) :: arr(:,:)

    type(c_ptr) :: ptr1, ptr2
    integer(c_size_t) :: ptr1_v, ptr2_v

    ! get full slices for all but the last slice
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2)+1))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)

    ret%ptr = ptr1_v
    ret%sz = (ptr2_v - ptr1_v) * (size(arr, 2) - 1)

    ! for the last slice, take the actual offset, otherwise we imght go OOB
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2)))
    ptr2 = c_loc(arr(lbound(arr,1)+1, lbound(arr,2)))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)
    ret%sz = ret%sz + (ptr2_v - ptr1_v) * size(arr, 1)
  end function
  function ext_acc_pass_3d_r4(arr) result(ret)
    implicit none
    type(ext_acc_arr_desc) :: ret
    real(4), intent(in) :: arr(:,:,:)

    type(c_ptr) :: ptr1, ptr2
    integer(c_size_t) :: ptr1_v, ptr2_v

    ! get full slices for all but the last slice
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr,3)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr,3)+1))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)

    ret%ptr = ptr1_v
    ret%sz = (ptr2_v - ptr1_v) * (size(arr, 3) - 1)

    ! for the last slice, take the actual offset, otherwise we imght go OOB
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr,3)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2)+1, lbound(arr,3)))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)
    ret%sz = ret%sz + (ptr2_v - ptr1_v) * size(arr, 2)
  end function
  function ext_acc_pass_4d_r4(arr) result(ret)
    implicit none
    type(ext_acc_arr_desc) :: ret
    real(4), intent(in) :: arr(:,:,:,:)

    type(c_ptr) :: ptr1, ptr2
    integer(c_size_t) :: ptr1_v, ptr2_v

    ! get full slices for all but the last slice
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3), lbound(arr,4)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3), lbound(arr,4)+1))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)

    ret%ptr = ptr1_v
    ret%sz = (ptr2_v - ptr1_v) * (size(arr, 4) - 1)

    ! for the last slice, take the actual offset, otherwise we imght go OOB
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3), lbound(arr,4)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3)+1, lbound(arr,4)))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)
    ret%sz = ret%sz + (ptr2_v - ptr1_v) * size(arr, 3)
  end function
  function ext_acc_pass_2d_r8(arr) result(ret)
    implicit none
    type(ext_acc_arr_desc) :: ret
    real(8), intent(in) :: arr(:,:)

    type(c_ptr) :: ptr1, ptr2
    integer(c_size_t) :: ptr1_v, ptr2_v

    ! get full slices for all but the last slice
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2)+1))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)

    ret%ptr = ptr1_v
    ret%sz = (ptr2_v - ptr1_v) * (size(arr, 2) - 1)

    ! for the last slice, take the actual offset, otherwise we imght go OOB
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2)))
    ptr2 = c_loc(arr(lbound(arr,1)+1, lbound(arr,2)))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)
    ret%sz = ret%sz + (ptr2_v - ptr1_v) * size(arr, 1)
  end function
  function ext_acc_pass_3d_r8(arr) result(ret)
    implicit none
    type(ext_acc_arr_desc) :: ret
    real(8), intent(in) :: arr(:,:,:)

    type(c_ptr) :: ptr1, ptr2
    integer(c_size_t) :: ptr1_v, ptr2_v

    ! get full slices for all but the last slice
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr,3)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr,3)+1))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)

    ret%ptr = ptr1_v
    ret%sz = (ptr2_v - ptr1_v) * (size(arr, 3) - 1)

    ! for the last slice, take the actual offset, otherwise we imght go OOB
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr,3)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2)+1, lbound(arr,3)))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)
    ret%sz = ret%sz + (ptr2_v - ptr1_v) * size(arr, 2)
  end function
  function ext_acc_pass_4d_r8(arr) result(ret)
    implicit none
    type(ext_acc_arr_desc) :: ret
    real(8), intent(in) :: arr(:,:,:,:)

    type(c_ptr) :: ptr1, ptr2
    integer(c_size_t) :: ptr1_v, ptr2_v

    ! get full slices for all but the last slice
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3), lbound(arr,4)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3), lbound(arr,4)+1))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)

    ret%ptr = ptr1_v
    ret%sz = (ptr2_v - ptr1_v) * (size(arr, 4) - 1)

    ! for the last slice, take the actual offset, otherwise we imght go OOB
    ptr1 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3), lbound(arr,4)))
    ptr2 = c_loc(arr(lbound(arr,1), lbound(arr,2), lbound(arr, 3)+1, lbound(arr,4)))
    ptr1_v= transfer(ptr1, ptr1_v)
    ptr2_v= transfer(ptr2, ptr2_v)
    ret%sz = ret%sz + (ptr2_v - ptr1_v) * size(arr, 3)
  end function
  function get_common_pointers(in_ptrs, out_ptrs) result(num_ranges)
    implicit none
    type(ext_acc_arr_desc), intent(in) :: in_ptrs(:)
    type(common_pointer_descr), intent(out) :: out_ptrs(:)

    integer(c_size_t), allocatable :: ptrs_only(:)
    logical, allocatable :: mask(:)
    integer, allocatable :: sort_index(:)

    type(ext_acc_arr_desc), allocatable :: common_ptrs(:)
    integer :: i, j, num_ranges
    integer(c_size_t) :: start1, start2, end1, end2
    logical :: found

    ! first sort the pointers increasingly such that no gaps are possible
    allocate(ptrs_only(size(in_ptrs)))
    do i = 1, size(in_ptrs)
      ptrs_only(i) = in_ptrs(i)%ptr
    enddo
    allocate(mask(size(in_ptrs)))
    do i = 1, size(in_ptrs)
      mask(i) = .true.
    enddo
    allocate(sort_index(size(in_ptrs)))
    do i = 1, size(in_ptrs)
      j = minloc(ptrs_only, 1, mask=mask)
      mask(j) = .false.
      sort_index(i) = j
    enddo

    ! initialize
    allocate(common_ptrs(size(in_ptrs)))
    do i = 1, size(in_ptrs)
      common_ptrs(1)%ptr = 0
      common_ptrs(1)%sz = 0
    enddo

    num_ranges = 1
    common_ptrs(1) = in_ptrs(sort_index(1))
    do i = 2, size(in_ptrs)
      found = .false.
      start1 = in_ptrs(sort_index(i))%ptr
      end1 = in_ptrs(sort_index(i))%ptr + in_ptrs(sort_index(i))%sz
      do j = 1, num_ranges
        start2 = common_ptrs(j)%ptr
        end2 = common_ptrs(j)%ptr + common_ptrs(j)%sz
        if (max(start1, start2) <= min(end1, end2)) then
          ! if we intersect with this range, extend the range
          common_ptrs(j)%ptr = min(start1, start2)
          common_ptrs(j)%sz = max(end1, end2) - common_ptrs(j)%ptr
          found = .true.
          exit
        endif
      enddo
      if (.not. found) then
        ! if we did not find anything: add a new one
        num_ranges = num_ranges + 1
        common_ptrs(num_ranges)%ptr = start1
        common_ptrs(num_ranges)%sz = end1 - start1
      endif
    enddo
    do i = 1, num_ranges
      out_ptrs(i)%ptr = transfer(common_ptrs(i)%ptr, out_ptrs(i)%ptr)
      out_ptrs(i)%sz = common_ptrs(i)%sz
    enddo
  end function
  subroutine ext_acc_create(ptrs, stream)
    use openacc
    implicit none
    type(ext_acc_arr_desc), intent(in) :: ptrs(:)
    integer(acc_handle_kind), optional :: stream

    type(common_pointer_descr), allocatable :: common_ptrs(:)

    integer :: i, num_ranges
    integer(4), pointer :: pp(:)
    integer(acc_handle_kind) :: stream_act

    if (present(stream)) then
      stream_act = stream
    else
      stream_act = acc_async_sync
    endif
    allocate(common_ptrs(size(ptrs)))
    num_ranges = get_common_pointers(ptrs, common_ptrs)

    do i = 1, num_ranges
      call c_f_pointer(common_ptrs(i)%ptr, pp, shape=[common_ptrs(i)%sz/sizeof(pp(1))])
      call acc_create_async(pp, common_ptrs(i)%sz, async=stream_act)
    enddo
  end subroutine
  subroutine ext_acc_copyin(ptrs, stream)
    use openacc
    implicit none
    type(ext_acc_arr_desc), intent(in) :: ptrs(:)
    integer(acc_handle_kind), optional :: stream

    type(common_pointer_descr), allocatable :: common_ptrs(:)

    integer :: i, num_ranges
    integer(4), pointer :: pp(:)

    integer(acc_handle_kind) :: stream_act

    if (present(stream)) then
      stream_act = stream
    else
      stream_act = acc_async_sync
    endif
    allocate(common_ptrs(size(ptrs)))
    num_ranges = get_common_pointers(ptrs, common_ptrs)

    do i = 1, num_ranges
      call c_f_pointer(common_ptrs(i)%ptr, pp, shape=[common_ptrs(i)%sz/sizeof(pp(1))])
      call acc_copyin_async(pp, common_ptrs(i)%sz, async=stream_act)
    enddo
  end subroutine
  subroutine ext_acc_copyout(ptrs, stream)
    use openacc
    implicit none
    type(ext_acc_arr_desc), intent(in) :: ptrs(:)
    integer(acc_handle_kind), optional :: stream

    type(common_pointer_descr), allocatable :: common_ptrs(:)

    integer :: i, num_ranges
    integer(4), pointer :: pp(:)

    integer(acc_handle_kind) :: stream_act

    if (present(stream)) then
      stream_act = stream
    else
      stream_act = acc_async_sync
    endif
    allocate(common_ptrs(size(ptrs)))
    num_ranges = get_common_pointers(ptrs, common_ptrs)

    do i = 1, num_ranges
      call c_f_pointer(common_ptrs(i)%ptr, pp, shape=[common_ptrs(i)%sz/sizeof(pp(1))])
      call acc_copyout_async(pp, common_ptrs(i)%sz, async=stream_act)
    enddo
  end subroutine
  subroutine ext_acc_delete(ptrs, stream)
    use openacc
    implicit none
    type(ext_acc_arr_desc), intent(in) :: ptrs(:)
    integer(acc_handle_kind), optional :: stream

    type(common_pointer_descr), allocatable :: common_ptrs(:)

    integer :: i, num_ranges
    integer(4), pointer :: pp(:)

    integer(acc_handle_kind) :: stream_act

    if (present(stream)) then
      stream_act = stream
    else
      stream_act = acc_async_sync
    endif
    allocate(common_ptrs(size(ptrs)))
    num_ranges = get_common_pointers(ptrs, common_ptrs)

    do i = 1, num_ranges
      call c_f_pointer(common_ptrs(i)%ptr, pp, shape=[common_ptrs(i)%sz/sizeof(pp(1))])
      call acc_delete_async(pp, common_ptrs(i)%sz, async=stream_act)
    enddo
  end subroutine
end module
