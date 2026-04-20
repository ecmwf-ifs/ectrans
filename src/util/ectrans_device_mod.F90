module ectrans_device_mod

implicit none
private

public :: ectrans_device_init, ectrans_device_is_host

contains

subroutine ectrans_device_init()
#if defined(_OPENACC)
  use openacc, only: acc_init, acc_get_device_type
  call acc_init(acc_get_device_type())
#endif
end subroutine

function ectrans_device_is_host() result(is_host)
#if defined(_OPENACC)
  use openacc, only: acc_get_device_type, acc_device_host
#elif defined(_OPENMP) && (_OPENMP >= 201307)
  ! OpenMP target device queries are part of the OpenMP 4.0 API.
  use omp_lib, only: omp_get_default_device, omp_get_initial_device
#endif
  logical :: is_host
  is_host = .true.
  call ectrans_device_init()
#if defined(_OPENACC)
  is_host = (acc_get_device_type() == acc_device_host)
#elif defined(_OPENMP) && (_OPENMP >= 201307)
  is_host = (omp_get_default_device() == omp_get_initial_device())
#endif
end function

!===================================================================================================

end module