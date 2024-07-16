! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

module device_mod

#ifdef CUDAGPU
#define hipDeviceSynchronize cudaDeviceSynchronize
#define hipStreamSynchronize cudaStreamSynchronize
#define hipStreamDestroy cudaStreamDestroy
#define hipSetDevice cudaSetDevice
#define hipGetDevice cudaGetDevice
#define hipGetDeviceCount cudaGetDeviceCount
#endif

interface device_sync

integer function device_synchronize() bind(C,name='hipDeviceSynchronize')
use iso_c_binding
end function device_synchronize

end interface device_sync

interface devicestreamsync

integer function device_stream_synchronize(stream) bind(C,name='hipStreamSynchronize')
use iso_c_binding
type(c_ptr) :: stream
end function device_stream_synchronize

end interface devicestreamsync

interface devicestreamdestroy

integer function device_stream_destroy(stream) bind(C,name='hipStreamDestroy')
use iso_c_binding
type(c_ptr) :: stream
end function device_stream_destroy

end interface devicestreamdestroy

interface devicesetdevice

integer function device_SetDevice(devnum) bind(C,name='hipSetDevice')
use iso_c_binding
integer(c_int),value :: devnum
end function device_SetDevice

end interface devicesetdevice

interface devicegetdevice

integer function device_GetDevice(devnum) bind(C,name='hipGetDevice')
use iso_c_binding
integer(c_int) :: devnum
end function device_GetDevice

end interface devicegetdevice

interface devicegetdevicecount

integer function device_GetDeviceCount(devnum) bind(C,name='hipGetDeviceCount')
use iso_c_binding
integer(c_int) :: devnum
end function device_GetDeviceCount

end interface devicegetdevicecount

interface devicegetmeminfo

integer function device_MemGetInfo(memfree_mb,memtotal_mb) bind(C,name='c_hipmemgetinfo')
use iso_c_binding
integer(c_int) :: memfree_mb, memtotal_mb
end function device_MemGetInfo

end interface devicegetmeminfo

end module device_mod
