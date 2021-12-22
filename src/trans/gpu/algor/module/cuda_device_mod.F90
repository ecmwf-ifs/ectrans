module cuda_device_mod

interface cuda_sync

integer function cuda_synchronize() bind(C,name='cudaDeviceSynchronize')
use iso_c_binding
end function cuda_synchronize

integer function cuda_stream_synchronize(stream) bind(C,name='cudaStreamSynchronize')
use iso_c_binding
type(c_ptr) :: stream
end function cuda_stream_synchronize

integer function cuda_stream_destroy(stream) bind(C,name='cudaStreamDestroy')
use iso_c_binding
type(c_ptr) :: stream
end function cuda_stream_destroy

end interface cuda_sync

interface cuda_device

integer function cuda_SetDevice(devnum) bind(C,name='cudaSetDevice')
use iso_c_binding
integer(c_int),value :: devnum
end function cuda_SetDevice

integer function cuda_GetDevice(devnum) bind(C,name='cudaGetDevice')
use iso_c_binding
integer(c_int) :: devnum
end function cuda_GetDevice

integer function cuda_GetDeviceCount(devnum) bind(C,name='cudaGetDeviceCount')
use iso_c_binding
integer(c_int) :: devnum
end function cuda_GetDeviceCount

end interface cuda_device

end module cuda_device_mod
