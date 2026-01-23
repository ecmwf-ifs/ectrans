module ectrans_memory
use, intrinsic :: iso_c_binding, only : c_char

private
public :: allocator

type allocator_t
contains
    procedure, nopass :: set_pinning
    procedure, nopass :: set_logging
    procedure, nopass :: set_logging_output_unit

    procedure, nopass, private :: allocate_var_real32_r1
    procedure, nopass, private :: allocate_var_real32_r2
    procedure, nopass, private :: allocate_var_real32_r3
    procedure, nopass, private :: allocate_var_real32_r4
    procedure, nopass, private :: allocate_var_real64_r1
    procedure, nopass, private :: allocate_var_real64_r2
    procedure, nopass, private :: allocate_var_real64_r3
    procedure, nopass, private :: allocate_var_real64_r4
    procedure, nopass, private :: allocate_var_label_real32_r1
    procedure, nopass, private :: allocate_var_label_real32_r2
    procedure, nopass, private :: allocate_var_label_real32_r3
    procedure, nopass, private :: allocate_var_label_real32_r4
    procedure, nopass, private :: allocate_var_label_real64_r1
    procedure, nopass, private :: allocate_var_label_real64_r2
    procedure, nopass, private :: allocate_var_label_real64_r3
    procedure, nopass, private :: allocate_var_label_real64_r4
    generic :: allocate => &
        & allocate_var_real32_r1, &
        & allocate_var_real32_r2, &
        & allocate_var_real32_r3, &
        & allocate_var_real32_r4, &
        & allocate_var_real64_r1, &
        & allocate_var_real64_r2, &
        & allocate_var_real64_r3, &
        & allocate_var_real64_r4, &
        & allocate_var_label_real32_r1, &
        & allocate_var_label_real32_r2, &
        & allocate_var_label_real32_r3, &
        & allocate_var_label_real32_r4, &
        & allocate_var_label_real64_r1, &
        & allocate_var_label_real64_r2, &
        & allocate_var_label_real64_r3, &
        & allocate_var_label_real64_r4

    procedure, nopass, private :: deallocate_var_real32_r1
    procedure, nopass, private :: deallocate_var_real32_r2
    procedure, nopass, private :: deallocate_var_real32_r3
    procedure, nopass, private :: deallocate_var_real32_r4
    procedure, nopass, private :: deallocate_var_real64_r1
    procedure, nopass, private :: deallocate_var_real64_r2
    procedure, nopass, private :: deallocate_var_real64_r3
    procedure, nopass, private :: deallocate_var_real64_r4
    procedure, nopass, private :: deallocate_var_label_real32_r1
    procedure, nopass, private :: deallocate_var_label_real32_r2
    procedure, nopass, private :: deallocate_var_label_real32_r3
    procedure, nopass, private :: deallocate_var_label_real32_r4
    procedure, nopass, private :: deallocate_var_label_real64_r1
    procedure, nopass, private :: deallocate_var_label_real64_r2
    procedure, nopass, private :: deallocate_var_label_real64_r3
    procedure, nopass, private :: deallocate_var_label_real64_r4
    generic :: deallocate => &
        & deallocate_var_real32_r1, &
        & deallocate_var_real32_r2, &
        & deallocate_var_real32_r3, &
        & deallocate_var_real32_r4, &
        & deallocate_var_real64_r1, &
        & deallocate_var_real64_r2, &
        & deallocate_var_real64_r3, &
        & deallocate_var_real64_r4, &
        & deallocate_var_label_real32_r1, &
        & deallocate_var_label_real32_r2, &
        & deallocate_var_label_real32_r3, &
        & deallocate_var_label_real32_r4, &
        & deallocate_var_label_real64_r1, &
        & deallocate_var_label_real64_r2, &
        & deallocate_var_label_real64_r3, &
        & deallocate_var_label_real64_r4
end type

type(allocator_t) :: allocator

character(kind=c_char), pointer, private :: c_label(:) => null()

interface
    function c_allocate_var(bytes) result(ptr) bind(c, name="ectrans_memory_allocate_var")
        use iso_c_binding, only: c_ptr, c_size_t
        integer(kind=c_size_t), value :: bytes
        type(c_ptr) :: ptr
    end function

    subroutine c_deallocate_var(ptr, bytes) bind(c, name="ectrans_memory_deallocate_var")
        use iso_c_binding, only: c_ptr, c_size_t
        type(c_ptr), value, intent(in) :: ptr
        integer(c_size_t), value, intent(in) :: bytes
    end subroutine

    subroutine c_set_pinning(pinning) bind(c, name="ectrans_memory_set_pinning")
        use iso_c_binding, only: c_int
        integer(c_int), value :: pinning
    end subroutine

    subroutine c_set_label(label) bind(c, name="ectrans_memory_set_label")
        use iso_c_binding, only: c_ptr
        type(c_ptr), value :: label ! must be null-terminated !
    end subroutine

    subroutine c_unset_label() bind(c, name="ectrans_memory_unset_label")
    end subroutine

    subroutine c_set_logging(logging) bind(c,name="ectrans_memory_set_logging")
        use iso_c_binding, only: c_int
        integer(c_int), value :: logging
    end subroutine

    subroutine c_set_logging_fortran_output_unit(output_unit) bind(c, &
          & name="ectrans_memory_set_logging_fortran_output_unit")
        use iso_c_binding, only: c_int
        integer(c_int), value :: output_unit
    end subroutine

end interface

contains

    function required_bytes(shape,real_kind)
        use iso_c_binding, only: c_int, c_size_t, c_float, c_double
        integer(c_int), intent(in) :: shape(:)
        integer, intent(in) :: real_kind
        integer(c_size_t) :: required_bytes
        if (real_kind == c_float) then
            required_bytes = product(int(shape,c_size_t)) * 4_c_size_t
        elseif (real_kind == c_double) then
            required_bytes = product(int(shape,c_size_t)) * 8_c_size_t
        else
            required_bytes = 0
        endif
    end function

    subroutine set_pinning(pinning)
        logical, intent(in) :: pinning
        if (pinning) then
            call c_set_pinning(1)
        else
            call c_set_pinning(0)
        endif
    end subroutine

    subroutine set_logging(logging)
        logical, intent(in) :: logging
        if (logging) then
            call c_set_logging(1)
        else
            call c_set_logging(0)
        endif
    end subroutine

    subroutine set_logging_output_unit(output_unit)
        integer, intent(in) :: output_unit
        call c_set_logging_fortran_output_unit(output_unit)
    end subroutine

    subroutine allocate_var_real32_r1(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        real(c_float), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_float))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine allocate_var_real32_r2(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        real(c_float), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_float))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine allocate_var_real32_r3(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_float))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine allocate_var_real32_r4(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_float))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine allocate_var_real64_r1(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        real(c_double), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_double))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine allocate_var_real64_r2(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        real(c_double), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_double))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine allocate_var_real64_r3(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_double))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine allocate_var_real64_r4(array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(:)
        type(c_ptr) :: mem
        mem = c_allocate_var(required_bytes(shape,c_double))
        call c_f_pointer(mem, array, shape)
    end subroutine

    subroutine deallocate_var_real32_r1(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc, c_null_ptr
        real(c_float), pointer, intent(inout) :: array(:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_float)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine deallocate_var_real32_r2(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc, c_null_ptr
        real(c_float), pointer, intent(inout) :: array(:,:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_float)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1,1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine deallocate_var_real32_r3(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc, c_null_ptr
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_float)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1,1,1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine deallocate_var_real32_r4(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc, c_null_ptr
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_float)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1,1,1,1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine deallocate_var_real64_r1(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc, c_null_ptr
        real(c_double), pointer, intent(inout) :: array(:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_double)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine deallocate_var_real64_r2(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc, c_null_ptr
        real(c_double), pointer, intent(inout) :: array(:,:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_double)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1,1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine deallocate_var_real64_r3(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc, c_null_ptr
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_double)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1,1,1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine deallocate_var_real64_r4(array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc, c_null_ptr
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        type(c_ptr) :: mem
        integer(c_size_t) :: bytes
        bytes = required_bytes(shape(array),c_double)
        mem = c_null_ptr
        if (bytes > 0) then
            mem = c_loc(array(1,1,1,1))
        endif
        call c_deallocate_var(mem, bytes)
        array => null()
    end subroutine

    subroutine set_label(label)
        use, intrinsic :: iso_c_binding, only : c_null_char, c_loc
        character(len=*), intent(in) :: label
        integer :: j, N
        if (associated(c_label)) then
            deallocate(c_label)
        endif
        N = len_trim(label)
        allocate(c_label(N+1))
        do j = 1, N
           c_label(j) = label(j:j)
        enddo
        c_label(N+1) = c_null_char
        call c_set_label(c_loc(c_label(1)))
    end subroutine

    subroutine unset_label()
        if (associated(c_label)) then
            deallocate(c_label)
        endif
        call c_unset_label()
    end subroutine

    subroutine allocate_var_label_real32_r1(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real32_r1(array, shape)
        call unset_label()
    end subroutine

    subroutine allocate_var_label_real32_r2(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real32_r2(array, shape)
        call unset_label()
    end subroutine

    subroutine allocate_var_label_real32_r3(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real32_r3(array, shape)
        call unset_label()
    end subroutine

    subroutine allocate_var_label_real32_r4(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real32_r4(array, shape)
        call unset_label()
    end subroutine

    subroutine allocate_var_label_real64_r1(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real64_r1(array, shape)
        call unset_label()
    end subroutine

    subroutine allocate_var_label_real64_r2(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real64_r2(array, shape)
        call unset_label()
    end subroutine

    subroutine allocate_var_label_real64_r3(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real64_r3(array, shape)
        call unset_label()
    end subroutine

    subroutine allocate_var_label_real64_r4(label, array, shape)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_int, c_f_pointer
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        integer(c_int), intent(in) :: shape(:)
        call set_label(label)
        call allocate_var_real64_r4(array, shape)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real32_r1(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:)
        call set_label(label)
        call deallocate_var_real32_r1(array)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real32_r2(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:)
        call set_label(label)
        call deallocate_var_real32_r2(array)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real32_r3(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:)
        call set_label(label)
        call deallocate_var_real32_r3(array)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real32_r4(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_float, c_loc
        character(len=*), intent(in) :: label
        real(c_float), pointer, intent(inout) :: array(:,:,:,:)
        call set_label(label)
        call deallocate_var_real32_r4(array)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real64_r1(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:)
        call set_label(label)
        call deallocate_var_real64_r1(array)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real64_r2(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:)
        call set_label(label)
        call deallocate_var_real64_r2(array)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real64_r3(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:)
        call set_label(label)
        call deallocate_var_real64_r3(array)
        call unset_label()
    end subroutine

    subroutine deallocate_var_label_real64_r4(label, array)
        use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_double, c_loc
        character(len=*), intent(in) :: label
        real(c_double), pointer, intent(inout) :: array(:,:,:,:)
        call set_label(label)
        call deallocate_var_real64_r4(array)
        call unset_label()
    end subroutine

    subroutine c_write_to_fortran_unit(unit,msg_cptr) bind(c, name="ectrans_memory_write_to_fortran_unit")
        use, intrinsic :: iso_c_binding, only: c_int32_t, c_ptr, c_char, c_associated
        integer(c_int32_t), value, intent(in) :: unit
        type(c_ptr), value, intent(in) :: msg_cptr
        character(kind=c_char,len=:), allocatable :: msg
        if( c_associated(msg_cptr) ) then
            call copy_c_ptr_to_string( msg_cptr, msg )
            write(unit,'(A)', advance='no') msg
        endif
    contains
        subroutine copy_c_str_to_string(s,string)
            use, intrinsic :: iso_c_binding
            character(kind=c_char,len=1), intent(in) :: s(:)
            character(len=:), allocatable :: string
            integer :: i, nchars
            do i = 1, size(s)
            if (s(i) == c_null_char) exit
            enddo
            nchars = i - 1  ! Exclude null character from Fortran string
            allocate( character(len=(nchars),kind=c_char) :: string )
            do i=1,nchars
            string(i:i) = s(i)
            enddo
        end subroutine
        subroutine copy_c_ptr_to_string(cptr,string)
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in) :: cptr
            character(kind=c_char,len=:), allocatable :: string
            character(kind=c_char), dimension(:), pointer  :: s
            integer(c_int), parameter :: MAX_STR_LEN = 2550
            call c_f_pointer ( cptr , s, (/MAX_STR_LEN/) )
            call copy_c_str_to_string( s, string )
        end subroutine
    end subroutine
end module
