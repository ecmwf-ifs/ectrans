! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program transform_test

!
! Spectral transform test
!
! This test performs spectral to real and real to spectral transforms repeated in
! timed loop.
!
!
! Author : George Mozdzynski
!

use parkind1, only: jpim, jprb, jprd
use oml_mod ,only : oml_max_threads
use mpl_mpif
use mpl_module
use yomgstats, only: jpmaxstat

implicit none

! Number of points in top/bottom latitudes
integer(kind=jpim), parameter :: min_octa_points = 20

integer(kind=jpim) :: istack, getstackusage
real(kind=jprb), dimension(1) :: zmaxerr(5), zerr(5)
real(kind=jprb) :: zmaxerrg

! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT
integer(kind=jpim), parameter :: noutdump = 7 ! Unit number for field output

! Default parameters
integer(kind=jpim) :: nlin   = 0   ! Linear grid (1) or not (0)
integer(kind=jpim) :: nq     = 2   ! Cubic grid (1) or cubic grid + collignon (2) or not (0)
integer(kind=jpim) :: nsmax  = 79  ! Spectral truncation
integer(kind=jpim) :: iters  = 10  ! Number of iterations for transform test
integer(kind=jpim) :: nflevg = 137 ! Default number of vertical levels

integer(kind=jpim) :: ndgl ! Number of latitudes
integer(kind=jpim) :: noutdum
integer(kind=jpim) :: nspec2
integer(kind=jpim) :: ngptot
integer(kind=jpim) :: ngptotg
integer(kind=jpim) :: ifld
integer(kind=jpim) :: iflds
integer(kind=jpim) :: icode
integer(kind=jpim) :: ioutsf
integer(kind=jpim) :: jroc
integer(kind=jpim) :: jb
integer(kind=jpim) :: ierr
integer(kind=jpim) :: nspec2g
integer(kind=jpim) :: iret
integer(kind=jpim) :: ntype
integer(kind=jpim) :: i
integer(kind=jpim) :: ja
integer(kind=jpim) :: ib
integer(kind=jpim) :: jprtrv

integer(kind=jpim), allocatable :: nloen(:), nprcids(:)
integer(kind=jpim) :: myproc, jj
integer :: jstep

real(kind=jprd) :: ztinit, ztloop, timef, ztstepmax, ztstepmin, ztstepavg, ztstepmed
real(kind=jprd) :: ztstepmax1, ztstepmin1, ztstepavg1, ztstepmed1
real(kind=jprd) :: ztstepmax2, ztstepmin2, ztstepavg2, ztstepmed2
real(kind=jprd), allocatable :: ztstep(:), ztstep1(:), ztstep2(:)

real(kind=jprb), allocatable :: znormsp(:), znormsp1(:), znormdiv(:), znormdiv1(:)
real(kind=jprb), allocatable :: znormvor(:), znormvor1(:), znormt(:), znormt1(:)
real(kind=jprb), allocatable :: znorm(:), znorm1(:)
real(kind=jprd) :: zaveave(0:jpmaxstat)

! Grid-point space data structures
real(kind=jprb), allocatable, target :: zgmv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), allocatable, target :: zgmvs  (:,:,:)   ! Single level fields at t and t-dt

! Spectral space data structures
real(kind=jprb), allocatable :: zspvorg(:,:)
real(kind=jprb), allocatable :: zspdivg(:,:)
real(kind=jprb), allocatable :: zspspg(:,:)
real(kind=jprb), allocatable :: zsptg(:,:,:)
real(kind=jprb), allocatable, target :: sp3d(:,:,:)
real(kind=jprb), pointer :: zvor(:,:) => null()
real(kind=jprb), pointer :: zdiv(:,:) => null()
real(kind=jprb), pointer :: zt(:,:,:) => null()
real(kind=jprb), allocatable :: zsp(:,:)

logical :: lstack = .false. ! Output stack info
logical :: luserpnm = .false.
logical :: lkeeprpnm = .false.
logical :: luseflt = .false. ! Use fast legendre transforms
logical :: ltrace_stats = .false.
logical :: lstats_omp = .false.
logical :: lstats_comms = .false.
logical :: lstats_mpl = .false.
logical :: lstats = .true. ! gstats statistics
logical :: lbarrier_stats = .false.
logical :: lbarrier_stats2 = .false.
logical :: ldetailed_stats = .false.
logical :: lstats_alloc = .false.
logical :: lsyncstats = .false.
logical :: lstatscpu = .false.
logical :: lstats_mem = .false.
logical :: lxml_stats = .false.
logical :: lfftw = .true. ! Use FFTW for Fourier transforms

integer(kind=jpim) :: nstats_mem = 0
integer(kind=jpim) :: ntrace_stats = 0
integer(kind=jpim) :: nprnt_stats = 1

logical :: lmpoff = .false. ! Message passing switch

! Whether to print verbose output or not
logical :: verbose = .false.

real(kind=jprb) :: zra=6371229._jprb

integer(kind=jpim) :: nmax_resol = 37 ! Max number of resolutions
integer(kind=jpim) :: npromatr = 0 ! nproma for trans lib
integer(kind=jpim) :: ncombflen = 1800000 ! Size of comm buffer

integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nthread
integer(kind=jpim) :: nprgpns ! Grid-point decomp
integer(kind=jpim) :: nprgpew ! Grid-point decomp
integer(kind=jpim) :: nprtrv ! Spectral decomp
integer(kind=jpim) :: nprtrw ! Spectral decomp
integer(kind=jpim) :: nspecresmin = 80 ! Minimum spectral resolution, for controlling nprtrw
integer(kind=jpim) :: mysetv
integer(kind=jpim) :: mysetw
integer(kind=jpim) :: mp_type = 2 ! Message passing type
integer(kind=jpim) :: mbx_size = 150000000 ! Mailbox size

integer(kind=jpim), allocatable :: numll(:), ivset(:)
integer(kind=jpim) :: ivsetsc(1)

integer(kind=jpim) :: nflevl

! sumpini
integer(kind=jpim) :: isqr
logical :: lsync_trans = .true. ! Activate barrier sync
logical :: leq_regions = .true. ! Eq regions flag


integer(kind=jpim) :: nproma
integer(kind=jpim) :: ngpblks
! locals
integer(kind=jpim) :: iprtrv
integer(kind=jpim) :: iprtrw
integer(kind=jpim) :: iprused, ilevpp, irest, ilev, jlev

logical :: llinfo

integer(kind=jpim) :: ndimgmv  = 9 ! Third dim. of gmv "(nproma,nflevg,ndimgmv,ngpblks)"
integer(kind=jpim) :: ndimgmvs = 3 ! Second dim. gmvs "(nproma,ndimgmvs,ngpblks)"

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "specnorm.h"
#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"

!===================================================================================================

call get_command_line_arguments(nsmax, iters, verbose)

write(nout,"(a)") "Running spectral transform test"
write(nout,"(a,i4)") "Resolution: TCO", nsmax
write(nout,"(a,i4)") "Iterations: ", iters

!===================================================================================================
! Message passing setup
! Participating processors limited by -P option
!===================================================================================================

call mpl_init()
!if( lstats ) call gstats(0,0)
ztinit=timef()

nproc= mpl_nproc()
myproc = mpl_myrank()
nthread= oml_max_threads()

! only output to stdout on pe 1
if( nproc > 1 ) then
  if( myproc /= 1 ) then
    open(unit=nout, file='/dev/null')
  endif
endif

if(ldetailed_stats)then
  lstats_omp=.true.
  lstats_comms=.true.
  lstats_mpl=.true.
  lstatscpu=.true.
  nprnt_stats=nproc
!  lstats_mem=.true.
!  lstats_alloc=.true.
endif

!===================================================================================================

allocate(nprcids(nproc))
do jj=1,nproc
  nprcids(jj) = jj
enddo

if( nproc <= 1 ) lmpoff=.true.
! Compute nprgpns and nprgpew
! This version selects most square-like distribution
! These will change if leq_regions=.true.
if( nproc == 0 ) nproc = 1
isqr=int(sqrt(real(nproc,jprb)))
do ja=isqr,nproc
  ib=nproc/ja
  if( ja*ib == nproc ) then
    nprgpns=max(ja,ib)
    nprgpew=min(ja,ib)
    exit
  endif
enddo

! From sumpini, although this
! should be specified in namelist
if( nspecresmin==0 ) nspecresmin=nproc

! Compute nprtrv and nprtrw
! if not provided in namelist
if( nprtrv > 0 .or. nprtrw > 0 ) then
  if( nprtrv == 0 ) nprtrv=nproc/nprtrw
  if( nprtrw == 0 ) nprtrw=nproc/nprtrv
  if( nprtrw*nprtrv /= nproc ) call abor1('transform_test:nprtrw*nprtrv /= nproc')
  if( nprtrw > nspecresmin ) call abor1('transform_test:nprtrw > nspecresmin')
else
  do jprtrv=4,nproc
    nprtrv=jprtrv
    nprtrw=nproc/nprtrv
    if( nprtrv*nprtrw /= nproc ) cycle
    if( nprtrv > nprtrw ) exit
    if( nprtrw > nspecresmin ) cycle
    if( nprtrw <= nspecresmin/(2*oml_max_threads()) ) exit
  enddo
  ! Go for approx square partition for backup
  if( nprtrv*nprtrw /= nproc .or. nprtrw > nspecresmin .or. nprtrv > nprtrw ) then
    isqr=int(sqrt(real(nproc,jprb)))
    do ja=isqr,nproc
      ib=nproc/ja
      if (ja*ib == nproc) then
        nprtrw=max(ja,ib)
        nprtrv=min(ja,ib)
        if (nprtrw > nspecresmin ) call abor1('transform_test:nprtrw &
                                           & (approx square value) > nspecresmin')
        exit
      endif
    enddo
  endif
endif

! Create communicators for mpi groups
if (.not.lmpoff) then
  call mpl_groups_create(nprtrw,nprtrv)
endif

if (lmpoff) then
  mysetw=(myproc-1)/nprtrv+1
  mysetv=mod(myproc-1,nprtrv)+1
else
  call mpl_cart_coords(myproc,mysetw,mysetv)
  ! Just checking for now...
  iprtrv=mod(myproc-1,nprtrv)+1
  iprtrw=(myproc-1)/nprtrv+1
  if (iprtrv/=mysetv .or. iprtrw/=mysetw) then
    call abor1('transform_test:inconsistency when computing mysetw and mysetv')
  endif
endif

if (.not.lmpoff) then
  llinfo=.false.
  if (myproc == 1) llinfo=.true.
  call mpl_buffer_method(kmp_type=mp_type,kmbx_size=mbx_size,kprocids=nprcids,ldinfo=llinfo)
endif

! Determine number of local levels for fourier and legendre calculations
! based on the values of nflevg and nprtrv
allocate(numll(nprtrv+1))

! Calculate remainder
iprused=min(nflevg+1,nprtrv)
ilevpp=nflevg/nprtrv
irest=nflevg-ilevpp*nprtrv
do jroc=1,nprtrv
  if(jroc <= irest) then
    numll(jroc)=ilevpp+1
  else
    numll(jroc)=ilevpp
  endif
enddo
numll(iprused+1:nprtrv+1)=0

nflevl=numll(mysetv)

ivsetsc(1)=iprused

ifld=0
iflds=0

icode = 0

!===================================================================================================
! Set resolution parameters
!===================================================================================================

! Calculate number of latitudes
ndgl = 2 * (nsmax + 1)

! Calculate number of points at each latitude for octahedral grid
allocate(nloen(ndgl))

do i = 1, ndgl / 2
  nloen(i) = min_octa_points + 4 * (i - 1)
  nloen(ndgl - i + 1) = nloen(i)
end do

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

call setup_trans0(kout=nout,kerr=nerr,kprintlev=merge(2, 0, verbose),kmax_resol=nmax_resol, &
&                 kpromatr=npromatr,kprgpns=nprgpns,kprgpew=nprgpew,kprtrw=nprtrw, &
&                 kcombflen=ncombflen,ldmpoff=lmpoff,ldsync_trans=lsync_trans, &
&                 ldeq_regions=leq_regions, &
&                 prad=zra,ldalloperm=.true.)

call setup_trans(ksmax=nsmax,kdgl=ndgl,kloen=nloen,ldsplit=.true.,&
&                 ldusefftw=.false., lduserpnm=luserpnm,ldkeeprpnm=lkeeprpnm, &
&                 lduseflt=luseflt)
!
call trans_inq(kspec2=nspec2,kspec2g=nspec2g,kgptot=ngptot,kgptotg=ngptotg)

! Default, no blocking
nproma=ngptot

! Calculate number of NPROMA blocks
ngpblks=(ngptot-1)/nproma+1

!===================================================================================================
! Initialize spectral arrays
!===================================================================================================

! Allocate spectral arrays
! Try to mimick IFS layout as much as possible
nullify(zvor)
nullify(zdiv)
nullify(zt)
allocate(sp3d(nflevl,nspec2,3))
allocate(zsp(1,nspec2))

call initialize_spectral_arrays(nsmax, zsp, sp3d)

! Point convenience variables to storage variable sp3d
zvor => sp3d(:,:,1)
zdiv => sp3d(:,:,2)
zt   => sp3d(:,:,3:3)

!===================================================================================================
! Print information before starting
!===================================================================================================

! Print configuration details
if (verbose .and. myproc == 1) then
  write(nout,'(a)')'===-=== Start of  runtime parameters ===-==='
  write(nout,'(" ")')
  write(nout,'("nlin=   ",i10)') nlin
  write(nout,'("nq=     ",i10)') nq
  write(nout,'("nsmax=  ",i10)') nsmax
  write(nout,'("ndgl=   ",i10)') ndgl
  write(nout,'("nproc=  ",i10)') nproc
  write(nout,'("nthread=",i10)') nthread
  write(nout,'("nprgpns=",i10)') nprgpns
  write(nout,'("nprgpew=",i10)') nprgpew
  write(nout,'("nprtrw= ",i10)') nprtrw
  write(nout,'("nprtrv= ",i10)') nprtrv
  write(nout,'("nproma= ",i10)') nproma
  write(nout,'("ngptot= ",i10)') ngptot
  write(nout,'("ngptotg=",i10)') ngptotg
  write(nout,'("nflevg= ",i10)') nflevg
  write(nout,'("iflds=  ",i10)') iflds
  write(nout,'("nspec2= ",i10)') nspec2
  write(nout,'("nspec2g=",i10)') nspec2g
  write(nout,'("luseflt=",l10)') luseflt
  write(nout,'(" ")')
  write(nout,'(a)') '===-=== End of   runtime parameters ===-==='
end if

allocate(ivset(nflevg))

! Compute spectral distribution
ilev = 0
do jb=1,nprtrv
  do jlev=1,numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo

! Allocate grid-point arrays
allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
allocate(zgmvs(nproma,ndimgmvs,ngpblks))

allocate(znormsp(1))
allocate(znormsp1(1))
allocate(znormvor(nflevg))
allocate(znormvor1(nflevg))
allocate(znormdiv(nflevg))
allocate(znormdiv1(nflevg))
allocate(znormt(nflevg))
allocate(znormt1(nflevg))

if( verbose ) then
  call specnorm(pspec=zvor(1:nflevl,:),pnorm=znormvor1,kvset=ivset(1:nflevg))
  call specnorm(pspec=zdiv(1:nflevl,:),pnorm=znormdiv1,kvset=ivset(1:nflevg))
  call specnorm(pspec=zt(1:nflevl,:,1),pnorm=znormt1,kvset=ivset(1:nflevg))
  call specnorm(pspec=zsp(1:1,:),      pnorm=znormsp1,kvset=ivsetsc(1:1))

  if(myproc == 1) then
    do ifld=1,1
      write(nout,'("sp  znorm(",i4,")=",f20.15)') ifld,znormsp1(ifld)
    enddo
    do ifld=1,nflevg
      write(nout,'("div znorm(",i4,")=",f20.15)') ifld,znormdiv1(ifld)
    enddo
    do ifld=1,nflevg
      write(nout,'("vor znorm(",i4,")=",f20.15)') ifld,znormvor1(ifld)
    enddo
    do ifld=1,nflevg
      write(nout,'("t   znorm(",i4,")=",f20.15)') ifld,znormt1(ifld)
    enddo
  endif
endif

ztinit=(timef()-ztinit)/1000.0_jprd
if (verbose .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a,i6,a,f9.2,a)') "transform_test initialisation, on",nproc,&
                                & " tasks, took",ztinit," sec"
  write(nout,'(" ")')
end if

if(iters<=0) call abor1('transform_test:iters <= 0')

allocate(ztstep(iters))
allocate(ztstep1(iters))
allocate(ztstep2(iters))

ztstepavg=0._jprd
ztstepmax=0._jprd
ztstepmin=9999999999999999._jprd
ztstepavg1=0._jprd
ztstepmax1=0._jprd
ztstepmin1=9999999999999999._jprd
ztstepavg2=0._jprd
ztstepmax2=0._jprd
ztstepmin2=9999999999999999._jprd

if (verbose .and. myproc == 1) then
  write(nout,'(a)') '===-=== start of  spec transforms  ===-==='
  write(nout,'(" ")')
end if

if( lstats ) then
  call gstats(0,0)
  call gstats_setup(nproc,myproc,nprcids,&
   & lstats,lstatscpu,lsyncstats,ldetailed_stats,lbarrier_stats,lbarrier_stats2,&
   & lstats_omp,lstats_comms,lstats_mem,nstats_mem,lstats_alloc,&
   & ltrace_stats,ntrace_stats,nprnt_stats,lxml_stats)
  call gstats_psut
  ! TODO: what is this?
  !call gstats_label_ifs
endif

ztloop=timef()

!===================================================================================================
! Do spectral transform loop
!===================================================================================================

do jstep=1,iters
  ztstep(jstep)=timef()

  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

  ztstep1(jstep)=timef()
  call inv_trans(pspvor=zvor,pspdiv=zdiv,pspsc2=zsp(1:1,:),&
     & pspsc3a=zt,&
     & ldscders=.true.,ldvorgp=.false.,lddivgp=.true.,lduvder=.false.,&
     & kresol=1,kproma=nproma,kvsetuv=ivset,kvsetsc2=ivsetsc(1:1),&
     & kvsetsc3a=ivset,&
     & pgpuv=zgmv(:,:,2:4,:),pgp2=zgmvs(:,1:3,:),&
     & pgp3a=zgmv(:,:,5:7,:))
  ztstep1(jstep)=(timef()-ztstep1(jstep))/1000.0_jprd

  !=================================================================================================
  ! While in grid point space, dump the values to disk
  !=================================================================================================

  ! dump a field to a binary file
  call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgmvs(:,1,:), 'S', noutdump)
  call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgmv(:,nflevg,3,:),  'U', noutdump)
  call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgmv(:,nflevg,4,:),  'V', noutdump)
  call dump_gridpoint_field(jstep, myproc, nproma, ngpblks, zgmv(:,nflevg,5,:),  'T', noutdump)

  !=================================================================================================
  ! Do direct transform
  !=================================================================================================

  ztstep2(jstep)=timef()
  call dir_trans(pspvor=zvor,pspdiv=zdiv,&
      & pspsc2=zsp(1:1,:),pspsc3a=zt,&
      & kresol=1,kproma=nproma,kvsetuv=ivset,kvsetsc2=ivsetsc(1:1),&
      & kvsetsc3a=ivset,&
      & pgpuv=zgmv(:,:,3:4,:),pgp2=zgmvs(:,1:1,:),&
      & pgp3a=zgmv(:,:,5:5,:))
  ztstep2(jstep)=(timef()-ztstep2(jstep))/1000.0_jprd

  !=================================================================================================
  ! Calculate timings
  !=================================================================================================

  ztstep(jstep)=(timef()-ztstep(jstep))/1000.0_jprd

  ztstepavg=ztstepavg+ztstep(jstep)
  ztstepmin=min(ztstep(jstep),ztstepmin)
  ztstepmax=max(ztstep(jstep),ztstepmax)

  ztstepavg1=ztstepavg1+ztstep1(jstep)
  ztstepmin1=min(ztstep1(jstep),ztstepmin1)
  ztstepmax1=max(ztstep1(jstep),ztstepmax1)

  ztstepavg2=ztstepavg2+ztstep2(jstep)
  ztstepmin2=min(ztstep2(jstep),ztstepmin2)
  ztstepmax2=max(ztstep2(jstep),ztstepmax2)

  !=================================================================================================
  ! Print norms
  !=================================================================================================

  if( verbose )then
    call specnorm(pspec=zsp(1:1,:),       pnorm=znormsp, kvset=ivsetsc(1:1))
    call specnorm(pspec=zvor(1:nflevl,:), pnorm=znormvor,kvset=ivset(1:nflevg))
    call specnorm(pspec=zdiv(1:nflevl,:), pnorm=znormdiv,kvset=ivset(1:nflevg))
    call specnorm(pspec=zt(1:nflevl,:,1), pnorm=znormt, kvset=ivset(1:nflevg))

    if( myproc==1 ) then
      ! Surface pressure
      zmaxerr(:)=-999.0
      do ifld=1,1
        zerr(1)=abs(znormsp1(ifld)/znormsp(ifld)-1.0_jprb)
        zmaxerr(1)=max(zmaxerr(1),zerr(1))
      enddo
      ! Divergence
      do ifld=1,nflevg
        zerr(2)=abs(znormdiv1(ifld)/znormdiv(ifld)-1.0_jprb)
        zmaxerr(2)=max(zmaxerr(2),zerr(2))
      enddo
      ! Vorticity
      do ifld=1,nflevg
        zerr(3)=abs(znormvor1(ifld)/znormvor(ifld)-1.0_jprb)
        zmaxerr(3)=max(zmaxerr(3),zerr(3))
      enddo
      ! Temperature
      do ifld=1,nflevg
        zerr(4)=abs(znormt1(ifld)/znormt(ifld)-1.0_jprb)
        zmaxerr(4)=max(zmaxerr(4),zerr(4))
      enddo
      write(nout,'("time step ",i6," took", f8.4," | sp max err="e10.3,&
                 & " | div max err="e10.3," | vor max err="e10.3," | t max err="e10.3)') &
                 &  jstep,ztstep(jstep),zmaxerr(1),zmaxerr(2),zmaxerr(3),zmaxerr(4)
    endif
  else
    if (verbose .and. myproc == 1) then
      write(nout,'("time step ",i6," took", f8.4)') jstep,ztstep(jstep)
    end if
  endif
enddo

ztloop=(timef()-ztloop)/1000.0_jprd

if (verbose .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a)') '===-=== End of   spec transforms  ===-==='
  write(nout,'(" ")')
end if


if( verbose ) then
  call specnorm(pspec=zvor(1:nflevl,:),pnorm=znormvor,kvset=ivset(1:nflevg))
  call specnorm(pspec=zdiv(1:nflevl,:),pnorm=znormdiv,kvset=ivset(1:nflevg))
  call specnorm(pspec=zt(1:nflevl,:,1),pnorm=znormt,kvset=ivset(1:nflevg))
  call specnorm(pspec=zsp(1:1,:),      pnorm=znormsp,kvset=ivsetsc(1:1))

  if(myproc == 1) then
    ! surface pressure
    zmaxerr(:)=-999.0
    do ifld=1,1
      zerr(1)=abs(znormsp1(ifld)/znormsp(ifld)-1.0d0)
      zmaxerr(1)=max(zmaxerr(1),zerr(1))
      write(nout,'("sp znorm(",i4,")=",f20.15," err=",e10.3)') ifld,znormsp(ifld),zerr(1)
    enddo
    ! divergence
    do ifld=1,nflevg
      zerr(2)=abs(znormdiv1(ifld)/znormdiv(ifld)-1.0d0)
      zmaxerr(2)=max(zmaxerr(2),zerr(2))
      write(nout,'("div znorm(",i4,")=",f20.15," err=",e10.3)') ifld,znormdiv(ifld),zerr(2)
    enddo
    ! vorticity
    do ifld=1,nflevg
      zerr(3)=abs(znormvor1(ifld)/znormvor(ifld)-1.0d0)
      zmaxerr(3)=max(zmaxerr(3),zerr(3))
      write(nout,'("vor znorm(",i4,")=",f20.15," err=",e10.3)') ifld,znormvor(ifld),zerr(3)
    enddo
    ! temperature
    do ifld=1,nflevg
      zerr(4)=abs(znormt1(ifld)/znormt(ifld)-1.0d0)
      zmaxerr(4)=max(zmaxerr(4),zerr(4))
      write(nout,'("t znorm(",i4,")=",f20.15," err=",e10.3)') ifld,znormt(ifld),zerr(4)
    enddo
    ! maximum error across all fields
    zmaxerrg=max(max(zmaxerr(1),zmaxerr(2)),max(zmaxerr(2),zmaxerr(3)))

    write(nout,'("surface pressure max error=",e10.3)')zmaxerr(1)
    write(nout,'("divergence       max error=",e10.3)')zmaxerr(2)
    write(nout,'("vorticity        max error=",e10.3)')zmaxerr(3)
    write(nout,'("temperature      max error=",e10.3)')zmaxerr(4)
    write(nout,'("global           max error=",e10.3)')zmaxerrg

  endif
endif

call mpl_allreduce(ztloop,     'sum', ldreprod=.false.)
call mpl_allreduce(ztstep,     'sum', ldreprod=.false.)
call mpl_allreduce(ztstepavg,  'sum', ldreprod=.false.)
call mpl_allreduce(ztstepmax,  'max', ldreprod=.false.)
call mpl_allreduce(ztstepmin,  'min', ldreprod=.false.)

call mpl_allreduce(ztstep1,    'sum', ldreprod=.false.)
call mpl_allreduce(ztstepavg1, 'sum', ldreprod=.false.)
call mpl_allreduce(ztstepmax1, 'max', ldreprod=.false.)
call mpl_allreduce(ztstepmin1, 'min', ldreprod=.false.)

call mpl_allreduce(ztstep2,    'sum', ldreprod=.false.)
call mpl_allreduce(ztstepavg2, 'sum', ldreprod=.false.)
call mpl_allreduce(ztstepmax2, 'max', ldreprod=.false.)
call mpl_allreduce(ztstepmin2, 'min', ldreprod=.false.)


ztstepavg=(ztstepavg/real(nproc,jprb))/real(iters,jprd)
ztloop=ztloop/real(nproc,jprd)
ztstep(:)=ztstep(:)/real(nproc,jprd)

call sort(ztstep,iters)
ztstepmed = ztstep(iters/2)

ztstepavg1=(ztstepavg1/real(nproc,jprb))/real(iters,jprd)
ztstep1(:)=ztstep1(:)/real(nproc,jprd)

call sort(ztstep1,iters)
ztstepmed1 = ztstep1(iters/2)

ztstepavg2=(ztstepavg2/real(nproc,jprb))/real(iters,jprd)
ztstep2(:)=ztstep2(:)/real(nproc,jprd)

call sort(ztstep2,iters)
ztstepmed2 = ztstep2(iters/2)

if (verbose .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a)') '===-=== Start   of time step stats ===-==='
  write(nout,'(" ")')
  write(nout,'("Inverse transforms")')
  write(nout,'("------------------")')
  write(nout,'("avg  (s): ",f8.4)') ztstepavg1
  write(nout,'("min  (s): ",f8.4)') ztstepmin1
  write(nout,'("max  (s): ",f8.4)') ztstepmax1
  write(nout,'("med  (s): ",f8.4)') ztstepmed1
  write(nout,'(" ")')
  write(nout,'("Direct transforms")')
  write(nout,'("-----------------")')
  write(nout,'("avg  (s): ",f8.4)') ztstepavg2
  write(nout,'("min  (s): ",f8.4)') ztstepmin2
  write(nout,'("max  (s): ",f8.4)') ztstepmax2
  write(nout,'("med  (s): ",f8.4)') ztstepmed2
  write(nout,'(" ")')
  write(nout,'("Inverse-direct transforms")')
  write(nout,'("-------------------------")')
  write(nout,'("avg  (s): ",f8.4)') ztstepavg
  write(nout,'("min  (s): ",f8.4)') ztstepmin
  write(nout,'("max  (s): ",f8.4)') ztstepmax
  write(nout,'("med  (s): ",f8.4)') ztstepmed
  write(nout,'("loop (s): ",f8.4)') ztloop
  write(nout,'(" ")')
  write(nout,'(a)') '===-=== End     of time step stats ===-==='
  write(nout,'(" ")')
endif

if( lstack ) then
!           Gather stack usage statistics
  istack = getstackusage()
  if(myproc == 1) then
    print 9000, istack
    9000 format("Stack utilisation information",/,&
         &"=============================",//,&
         &"Task           size(bytes)",/,&
         &"====           ===========",//,&
         &"   1",11x,i10)

    do i=2,nproc
      call mpl_recv(istack,ksource=nprcids(i),ktag=i, &
           & cdstring='transform_test:')
      print '(i4,11x,i10)', i,istack
    enddo
  else
    call mpl_send(istack,kdest=nprcids(1),ktag=myproc, &
             &   cdstring='transform_test:')
  endif
endif


!===================================================================================================

if( lstats ) then
  call gstats(0,1)
  call gstats_print(nout,zaveave,jpmaxstat)
endif

!===================================================================================================

! Close file
if( nproc > 1 ) then
  if( myproc /= 1 ) then
    close(unit=nout)
  endif
endif

deallocate(zgmv)
deallocate(zgmvs)

!===================================================================================================
! Finalize MPI
!===================================================================================================

call mpl_barrier()
call mpl_end()

!===================================================================================================

contains

subroutine get_command_line_arguments(nsmax, iters, verbose)

  integer, intent(inout) :: nsmax   ! Spectral truncation
  integer, intent(inout) :: iters   ! Number of iterations for transform test
  logical, intent(inout) :: verbose ! Print verbose output or not

  character(len=128) :: carg          ! Storage variable for command line arguments
  integer            :: iarg = 1      ! Argument index
  integer            :: verbosity = 0 ! Level of verbosity (0, 1 or 2)
  integer            :: stat          ! For storing success status of string->integer conversion

  do while (iarg <= command_argument_count())
    call get_command_argument(iarg, carg)

    select case(carg)
      ! Parse help argument
      case('-h')
        call print_help()
        stop
      ! Parse verbosity argument
      case('-v')
        verbosity = 1
      ! Parse number of iterations argument
      case('-n')
        iarg = iarg + 1
        call get_command_argument(iarg, carg)
        call str2int(carg, iters, stat)
        if (stat /= 0 .or. iters < 1) then
          call print_help()
          call abor1("Invalid argument for -n: " // carg)
        end if
      ! Parse spectral truncation argument
      case('-t')
        iarg = iarg + 1
        call get_command_argument(iarg, carg)
        call str2int(carg, nsmax, stat)
        if (stat /= 0 .or. nsmax < 1) then
          call print_help()
          call abor1("Invalid argument for -t: " // carg)
        end if
    end select
    iarg = iarg + 1
  end do

  ! TODO: implement different levels of verbosity
  verbose = verbosity == 1

end subroutine get_command_line_arguments

!===================================================================================================

subroutine str2int(str, int, stat)
  character(len=*), intent(in) :: str
  integer, intent(out) :: int
  integer, intent(out) :: stat

  read(str, *, iostat=stat) int
end subroutine str2int

!===================================================================================================

subroutine sort(a, n)

  real(kind=jprd), intent(inout) :: a(n)
  integer(kind=jpim), intent(in) :: n

  real(kind=jprd) :: x

  integer :: i, j

  do i = 2, n
    x = a(i)
    j = i - 1
    do while (j >= 1)
      if (a(j) <= x) exit
        a(j + 1) = a(j)
        j = j - 1
      end do
    a(j + 1) = x
  end do

end subroutine sort

!===================================================================================================

subroutine print_help

  write(nout, "(a)") ""

  if (jprb == jprd) then
    write(nout, "(a)") "NAME    driver-spectrans-dp"
  else
    write(nout, "(a)") "NAME    driver-spectrans-sp"
  end if

  write(nout, "(a)") "DESCRIPTION"
  write(nout, "(a)") "        This program tests ecTrans by transforming fields back and forth&
    &between spectral "
  if (jprb == jprd) then
    write(nout, "(a)") "        space and grid-point space (double-precision version)"
  else
    write(nout, "(a)") "        space and grid-point space (single-precision version)"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "USAGE"
  if (jprb == jprd) then
    write(nout, "(a)") "        driver-spectrans-dp [options]"
  else
    write(nout, "(a)") "        driver-spectrans-sp [options]"
  end if

  write(nout, "(a)") "OPTIONS"
  write(nout, "(a)") "        -h      Print this message"
  write(nout, "(a)") "        -n iterations"
  write(nout, "(a)") "                Run for this many inverse/direct transform iterations"
  write(nout, "(a)") "        -t truncation"
  write(nout, "(a)") "                Run with this triangular spectral truncation"

  write(nout, "(a)") ""
  write(nout, "(a)") ""

end subroutine print_help

!===================================================================================================

subroutine initialize_spectral_arrays(nsmax, zsp, sp3d)

  integer,         intent(in)    :: nsmax       ! Spectral truncation
  real(kind=jprb), intent(inout) :: zsp(:,:)    ! Surface pressure
  real(kind=jprb), intent(inout) :: sp3d(:,:,:) ! 3D fields

  integer(kind=jpim) :: nflevl
  integer(kind=jpim) :: nfield

  integer :: i, j

  nflevl = size(sp3d, 1)
  nfield = size(sp3d, 3)

  ! First initialize surface pressure
  call initialize_2d_spectral_field(nsmax, zsp(1,:))

  ! Then initialize all of the 3D fields
  do i = 1, nflevl
    do j = 1, nfield
      call initialize_2d_spectral_field(nsmax, sp3d(i,:,j))
    end do
  end do

end subroutine initialize_spectral_arrays

!===================================================================================================

subroutine initialize_2d_spectral_field(nsmax, field)

  integer,         intent(in)    :: nsmax    ! Spectral truncation
  real(kind=jprb), intent(inout) :: field(:) ! Field to initialize

  integer :: i, index, num_my_zon_wns
  integer, allocatable :: my_zon_wns(:), nasm0(:)

  ! Choose a spherical harmonic to initialize arrays
  integer :: m_num = 4  ! Zonal wavenumber
  integer :: l_num = 19  ! Total wavenumber

  ! First initialise all spectral coefficients to zero
  field(:) = 0.0

  ! Get zonal wavenumbers this rank is responsible for
  call trans_inq(knump=num_my_zon_wns)
  allocate(my_zon_wns(num_my_zon_wns))
  call trans_inq(kmyms=my_zon_wns)

  ! If rank is responsible for the chosen zonal wavenumber...
  if (any(my_zon_wns == m_num) ) then
    ! Get array of spectral array addresses (this maps (m, n=m) to array index)
    allocate(nasm0(0:nsmax))
    call trans_inq(kasm0=nasm0)

    ! Find out local array index of chosen spherical harmonic
    index = nasm0(m_num) + 2 * (l_num - m_num) + 1

    ! Set just that element to a constant value
    field(index) = 1.0
  else
    return
  end if

end subroutine initialize_2d_spectral_field

!===================================================================================================

subroutine dump_gridpoint_field(jstep, myproc, nproma, ngpblks, fld, fldchar, noutdump)

! Dump a 2d field to a binary file.

integer(kind=jpim), intent(in) :: jstep ! Time step, used for naming file
integer(kind=jpim), intent(in) :: myproc ! MPI rank, used for naming file
integer(kind=jpim), intent(in) :: nproma ! Size of nproma
integer(kind=jpim), intent(in) :: ngpblks ! Number of nproma blocks
real(kind=jprb)   , intent(in) :: fld(nproma,ngpblks) ! 2D field
character         , intent(in) :: fldchar ! Single character field identifier
integer(kind=jpim), intent(in) :: noutdump ! Tnit number for output file

character(len=14) :: filename = "x.xxx.xxxx.dat"

write(filename(1:1),'(a1)') fldchar
write(filename(3:5),'(i3.3)') jstep
write(filename(7:10),'(i4.4)') myproc

open(noutdump, file=filename, form="unformatted")
write(noutdump) reshape(fld, (/ nproma*ngpblks /))
close(noutdump)

end subroutine dump_gridpoint_field

end program transform_test

!===================================================================================================
