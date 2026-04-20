module ectrans_grids_mod
implicit none
private
public :: cubic_octahedral_gaussian_grid
public :: parse_gaussian_grid

contains

subroutine parse_gaussian_grid(cgrid,ndgl,nloen)
  use iso_fortran_env, only : int32
  use ectrans_error_mod, only : ectrans_error_parsing_failed
  implicit none
  character(len=*), intent(in) :: cgrid
  integer(kind=int32), intent(inout) :: ndgl
  integer(kind=int32), intent(inout), allocatable :: nloen(:)

  integer :: ios
  integer(kind=int32) :: gaussian_number
  integer(kind=int32) :: i

  read(cgrid(2:len_trim(cgrid)),*,IOSTAT=ios) gaussian_number
  if (ios==0) then
    ndgl = 2 * gaussian_number
    allocate(nloen(ndgl))
    if (cgrid(1:1) == 'F') then ! Regular Gaussian grid
      nloen(:) = gaussian_number * 4
      return
    endif
    if (cgrid(1:1) == 'O') then ! Octahedral Gaussian grid
      do i = 1, ndgl / 2
        nloen(i) = 20 + 4 * (i - 1)
        nloen(ndgl - i + 1) = nloen(i)
      end do
      return
    endif
  endif
  call ectrans_error_parsing_failed("ERROR: Unsupported grid specified: "// trim(cgrid), print_help)
contains
  subroutine print_help(unit)
    integer, intent(in) :: unit
    write(unit,'(A)') "Supported grids: F<N>, O<N> where N is the number of Gaussian latitudes between pole and equator" // &
      & " (e.g., F128, O128)"
  end subroutine
end subroutine

function cubic_octahedral_gaussian_grid(nsmax) result(cgrid)
  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'O',nsmax+1
end function

end module