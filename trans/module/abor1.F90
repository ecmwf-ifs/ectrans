subroutine abor1(cl)
use tpm_gen
use tpm_distr
implicit none
character*(*) cl
integer ierror
close(nout)
write(0,'(a)') cl
if(nproc > 1 ) then
 call mpe_barrier(ierror)
 call mpe_end(ierror)
endif
stop
end
