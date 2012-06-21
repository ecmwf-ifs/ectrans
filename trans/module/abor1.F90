subroutine abor1(cl)
USE MPL_MODULE
use tpm_gen
use tpm_distr
implicit none
character*(*) cl
integer ierror
close(nout)
write(0,'(a)') cl
if(nproc > 1 ) then
 call mpl_barrier()
 call mpl_end()
endif
stop
end
