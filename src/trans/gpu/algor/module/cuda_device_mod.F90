! (C) Copyright 2020- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

module cuda_device_mod

interface device_sync

integer function device_synchronize() bind(C,name='cudaDeviceSynchronize')
use iso_c_binding
end function device_synchronize

end interface device_sync

interface devicestreamsync

integer function device_stream_synchronize(stream) bind(C,name='cudaStreamSynchronize')
use iso_c_binding
type(c_ptr) :: stream
end function device_stream_synchronize

end interface devicestreamsync

interface devicestreamdestroy

integer function device_stream_destroy(stream) bind(C,name='cudaStreamDestroy')
use iso_c_binding
type(c_ptr) :: stream
end function device_stream_destroy

end interface devicestreamdestroy

interface devicesetdevice

integer function device_SetDevice(devnum) bind(C,name='cudaSetDevice')
use iso_c_binding
integer(c_int),value :: devnum
end function device_SetDevice

end interface devicesetdevice

interface devicegetdevice

integer function device_GetDevice(devnum) bind(C,name='cudaGetDevice')
use iso_c_binding
integer(c_int) :: devnum
end function device_GetDevice

end interface devicegetdevice

interface devicegetdevicecount

integer function device_GetDeviceCount(devnum) bind(C,name='cudaGetDeviceCount')
use iso_c_binding
integer(c_int) :: devnum
end function device_GetDeviceCount

end interface devicegetdevicecount

end module cuda_device_mod
