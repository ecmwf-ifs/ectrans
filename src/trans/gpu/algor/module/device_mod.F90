! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

module device_mod

#ifdef HIPGPU
use hip_device_mod
#endif
#ifdef CUDAGPU
use cuda_device_mod
#endif

implicit none

end module device_mod
