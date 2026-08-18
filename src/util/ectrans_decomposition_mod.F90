module ectrans_decomposition_mod
use iso_fortran_env, only : int32, real64
implicit none
private

public :: ectrans_spectral_decomposition
public :: ectrans_gridpoint_decomposition
public :: ectrans_make_spectral_distribution

contains

subroutine get_square_like_gridpoint_decomposition(nproc, kprgpew, kprgpns)
  use iso_fortran_env, only : real64
  ! Try to find a good 2D distribution of processes over the grid points
  ! This version selects most square-like distribution
  ! These will change if leq_regiontrass=.true.
  integer(kind=int32), intent(in) :: nproc
  integer(kind=int32), intent(out) :: kprgpew, kprgpns
  integer(kind=int32) :: isqr, ja, ib
  isqr = int(sqrt(real(nproc,real64)))
  do ja = isqr, nproc
    ib = nproc/ja
    if (ja*ib == nproc) then
      kprgpns = max(ja,ib)
      kprgpew = min(ja,ib)
      exit
    endif
  enddo
end subroutine

subroutine ectrans_gridpoint_decomposition(nproc, kprgpns, kprgpew)
  integer(kind=int32), intent(in) :: nproc
  integer(kind=int32), intent(inout) :: kprgpns, kprgpew
  if (nproc == 1) then
    kprgpns = 1
    kprgpew = 1
    return
  endif
  call get_square_like_gridpoint_decomposition(nproc, kprgpew, kprgpns)
end subroutine

subroutine ectrans_spectral_decomposition(nproc, kprtrw, kprtrv)
  integer(kind=int32), intent(in) :: nproc
  integer(kind=int32), intent(inout) :: kprtrw, kprtrv
  integer(kind=int32) :: jprtrv, isqr, ja, ib, nprtrv, nprtrw
  if (nproc == 1) then
    kprtrw = 1
    kprtrv = 1
    return
  endif
  nprtrw = kprtrw
  nprtrv = kprtrv
  if (nprtrw > 0 .or. nprtrv > 0) then
    if (nprtrw == 0) nprtrw = nproc/nprtrv
    if (nprtrv == 0) nprtrv = nproc/nprtrw
    if (nprtrw*nprtrv /= nproc) call abor1('nprtrw*nprtrv /= nproc')
  else
    do jprtrv = 4, nproc
      nprtrv = jprtrv
      nprtrw = nproc/nprtrv
      if (nprtrv*nprtrw /= nproc) cycle
      if (nprtrv > nprtrw) exit
    enddo
    ! Go for approx square partition for backup
    if (nprtrv*nprtrw /= nproc .or. nprtrv > nprtrw) then
      isqr = int(sqrt(real(nproc,real64)),int32)
      do ja = isqr, nproc
        ib = nproc/ja
        if (ja*ib == nproc) then
          nprtrw = max(ja, ib)
          nprtrv = min(ja, ib)
          exit
        endif
      enddo
    endif
  endif
  kprtrw = nprtrw
  kprtrv = nprtrv
end subroutine

function ectrans_make_spectral_distribution(nfld, kprtrv) result(num_fields_per_proc)
  use ec_parkind, only : jpim
  implicit none
  integer, intent(in) :: nfld
  integer(kind=jpim), intent(in) :: kprtrv
  integer(kind=jpim), allocatable :: num_fields_per_proc(:)
  integer(kind=jpim) :: ilevpp, irest, jroc
  allocate(num_fields_per_proc(kprtrv))
  ilevpp = nfld/kprtrv
  irest = nfld - ilevpp*kprtrv
  do jroc = 1, kprtrv
    if (jroc <= irest) then
      num_fields_per_proc(jroc) = ilevpp+1
    else
      num_fields_per_proc(jroc) = ilevpp
    endif
  enddo
end function

end module