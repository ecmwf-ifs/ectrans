module sutrle_mod
contains
subroutine sutrle(pnm)

!**** *sutrle * - transposition of Legendre polynomials during set-up

!     Purpose.
!     --------
!           transposition of Legendre polynomials during set-up

!**   Interface.
!     ----------
!        *call* *sutrle(pnm)

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!     ------------------------------------------------------------------


#include "tsmbkind.h"

use tpm_gen
use tpm_dim
use tpm_distr
use tpm_fields
use set2pe_mod


IMPLICIT NONE

REAL_B :: pnm(d%nlei3d,r%nspoleg)

REAL_B, allocatable :: zcombuf(:)
REAL_B, pointer     :: zpnm(:,:)
!     LOCAL INTEGER SCALARS
INTEGER_M :: ierr, iglloc, ilrec, im, inentr, ipos, irecd,&
             &irecset, irecv, isend, isendset, itag, itagr, &
             &jgl, jglloc, jm, jmloc, jn, jroc ,iofft, ioffg

!     LOCAL LOGICAL SCALARS
LOGICAL :: lladmsg, llexact

!     ------------------------------------------------------------------

!*       0.    Some initializations.
!              ---------------------

!! Workaround for obscure unwillingness to vectorize on VPP
zpnm => f%rpnm

! Perform barrier synchronisation to guarantee all processors have
! completed all previous communication

if( nproc > 1 )then
  call mpe_barrier(ierr)
  if( ierr /= 0 )then
    call abor1('sutrle: at start, error in mpe_barrier')
  endif
endif

allocate (zcombuf(ncombflen))

do jroc=1,nprtrw-1

  lladmsg = .false.
  itag = mtagletr

!*     Define PE to which data have to be sent and PE from which
!*     data have to be received

  isend = mysetw-jroc
  irecv = mysetw+jroc
  if (isend <= 0)     isend = isend+nprtrw
  if (irecv > nprtrw) irecv = irecv-nprtrw
  irecset = irecv
  isendset = isend
  call set2pe(isend,0,0,isend,mysetv)
  call set2pe(irecv,0,0,irecv,mysetv)

!*   copy data to be sent into zcombuf

  ipos=0
  do jm=0,r%nsmax
    if (isendset == d%nprocm(jm)) then
      inentr = (d%nlatle(mysetw)-d%nlatls(mysetw)+1)*(r%ntmax-jm+2)
      if (ipos + inentr < ncombflen) then
        do jgl=d%nlatls(mysetw),d%nlatle(mysetw)
          jglloc = jgl - d%nlatls(mysetw) + 1
          do jn=1,r%ntmax-jm+2
            ipos = ipos + 1
            zcombuf(ipos)=pnm(jglloc,d%npmg(jm)+jn)
          enddo
        enddo
      else
        do jgl=d%nlatls(mysetw),d%nlatle(mysetw)
          jglloc = jgl - d%nlatls(mysetw) + 1
          do jn=1,r%ntmax-jm+2
            ipos = ipos + 1
            zcombuf(ipos)=pnm(jglloc,d%npmg(jm)+jn)
            if (ipos == ncombflen) then
              call mpe_send(zcombuf,ipos,mrealt,&
               &nprcids(isend),itag,0,0,0,ierr)
              ipos = 0
              itag = itag + 1
              llexact=jgl == d%nlatle(mysetw).and.jn == r%ntmax-jm+2
              if (.not.llexact) lladmsg = .true.
            endif
          enddo
        enddo
      endif
    endif
  enddo

!*   send message (if not empty or if message has been split)

  if (ipos > 0.or. lladmsg) then
    call mpe_send(zcombuf,ipos,mrealt,nprcids(isend),itag,0,0,0,ierr)
  endif

  ilrec = 0
  itag = mtagletr
  if (d%nump > 0.and. d%nlatle(irecset) >= d%nlatls(irecset)) then

!*   receive message (if not empty)

    call mpe_recv(zcombuf,ncombflen,mrealt,nprcids(irecv),itag,&
     &0,0,0,ilrec,irecd,itagr,ierr)

!*   copy data from buffer to f%rpnm

    ipos=0
    do jmloc=1,d%nump
      jm = d%myms(jmloc)
      inentr = (d%nlatle(irecset)-d%nlatls(irecset)+1)*(r%ntmax-jm+2)
      iofft = d%npmt(jm) 
      if (ipos + inentr < ncombflen) then
        do jgl=d%nlatls(irecset),d%nlatle(irecset)
          do jn=1,r%ntmax-jm+2
            ipos = ipos + 1
            zpnm(jgl,iofft+jn) = zcombuf(ipos)
          enddo
        enddo
      else
        do jgl=d%nlatls(irecset),d%nlatle(irecset)
          do jn=1,r%ntmax-jm+2
            ipos = ipos + 1
            zpnm(jgl,iofft+jn) = zcombuf(ipos)
            if (ipos == ncombflen) then
              itag = itag + 1
              call mpe_recv(zcombuf,ncombflen,mrealt,&
               &nprcids(irecv),itag,0,0,0,ilrec,irecd,itagr,&
               &ierr)
              ipos = 0
            endif
          enddo
        enddo
      endif
    enddo

!*    check received message length

    if (ilrec /= ipos) then
      write(nout,*)' sutrle: ilrec,ipos,ncomblen ',ilrec,ipos,ncombflen
      call abor1(' sutrle:received message length does not match')
    endif
  endif

! Perform barrier synchronisation to guarantee all processors have
! completed communication for this jroc loop iteration

  call mpe_barrier(ierr)
  if( ierr /= 0 )then
    call abor1('sutrle: jroc loop, error in mpe_barrier')
  endif

enddo

!*    copy data from pnm to rpnm

!$OMP PARALLEL PRIVATE(jmloc,im,jgl,iglloc,jn)
!$OMP DO SCHEDULE(STATIC,1)
do jmloc=1,d%nump
  im = d%myms(jmloc)
  iofft = d%npmt(im)
  ioffg = d%npmg(im)
  do jgl=d%nlatls(mysetw),d%nlatle(mysetw)
    iglloc = jgl-d%nlatls(mysetw)+1
    do jn=1,r%ntmax-im+2
      zpnm(jgl,iofft+jn) = pnm(iglloc,ioffg+jn)
    enddo
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

deallocate (zcombuf)

end subroutine sutrle
end module sutrle_mod
