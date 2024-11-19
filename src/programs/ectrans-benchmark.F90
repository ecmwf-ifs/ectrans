! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program ectrans_benchmark

#ifdef USE_PINNED
#define PINNED_TAG , pinned
#else
#define PINNED_TAG
#endif

!
! Spectral transform test
!
! This test performs spectral to real and real to spectral transforms repeated in
! timed loop.
!
! Authors : George Mozdzynski
!           Willem Deconinck
!           Ioan Hadade
!           Sam Hatfield
!

use parkind1, only: jpim, jprb, jprd
use oml_mod ,only : oml_max_threads
use mpl_module
use yomgstats, only: jpmaxstat, gstats_lstats => lstats
use yomhook, only : dr_hook_init

implicit none

! Number of points in top/bottom latitudes
integer(kind=jpim), parameter :: min_octa_points = 20

integer(kind=jpim) :: istack, getstackusage
real(kind=jprd) :: zmaxerr(5)
real(kind=jprd) :: zmaxerrg

! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT
integer(kind=jpim), parameter :: noutdump = 7 ! Unit number for field output

! Default parameters
integer(kind=jpim) :: iters   = 10  ! Number of iterations for transform test
integer(kind=jpim) :: nfld    = 1   ! Number of 3D scalar fields
integer(kind=jpim) :: nlev    = 1   ! Number of vertical levels
integer(kind=jpim) :: iters_warmup = 3 ! Number of warm up steps (for which timing statistics should be ignored)

integer(kind=jpim) :: nflevg  ! Total number of vertical levels

integer(kind=jpim) :: nspec2  ! Number of spectral coefficients (real and imaginary)
integer(kind=jpim) :: ngptot  ! Total number of grid points on this task
integer(kind=jpim) :: ngptotg ! Total number of grid points across all tasks

integer(kind=jpim) :: ifld
integer(kind=jpim) :: jroc
integer(kind=jpim) :: jb
integer(kind=jpim) :: nspec2g
integer(kind=jpim) :: i
integer(kind=jpim) :: ja
integer(kind=jpim) :: ib
integer(kind=jpim) :: jprtrv

integer(kind=jpim), allocatable :: nprcids(:)
integer(kind=jpim) :: myproc, jj
integer :: jstep

real(kind=jprd) :: ztinit, ztloop, timef, ztstepmax, ztstepmin, ztstepavg, ztstepmed
real(kind=jprd) :: ztstepmax1, ztstepmin1, ztstepavg1, ztstepmed1
real(kind=jprd) :: ztstepmax2, ztstepmin2, ztstepavg2, ztstepmed2
real(kind=jprd), allocatable :: ztstep(:), ztstep1(:), ztstep2(:)

real(kind=jprb), allocatable :: znormvor(:), znormvor1(:), znormdiv(:), znormdiv1(:)
real(kind=jprb), allocatable :: znormscalar(:), znormscalar1(:)
real(kind=jprb), allocatable :: znormsc3a(:), znormsc3a1(:), znormsc2(:), znormsc21(:)

real(kind=jprd) :: zaveave(0:jpmaxstat)

! Spectral space data structures
real(kind=jprb), allocatable PINNED_TAG :: zspvor(:,:)
real(kind=jprb), allocatable PINNED_TAG :: zspdiv(:,:)
real(kind=jprb), allocatable PINNED_TAG :: zspscalar(:,:)
real(kind=jprb), allocatable PINNED_TAG :: zspsc3a(:,:,:)
real(kind=jprb), allocatable PINNED_TAG :: zspsc2(:,:)

! Grid-point space data structures
real(kind=jprb), allocatable PINNED_TAG :: zgp(:,:,:)
real(kind=jprb), allocatable PINNED_TAG :: zgpuv(:,:,:,:)
real(kind=jprb), allocatable PINNED_TAG :: zgp3a(:,:,:,:)
real(kind=jprb), allocatable PINNED_TAG :: zgp2(:,:,:)

logical :: lstack = .false. ! Output stack info

! setup_trans options
integer(kind=jpim) :: nsmax   = 79  ! Spectral truncation
integer(kind=jpim) :: ndgl    ! Number of latitudes
integer(kind=jpim), allocatable :: nloen(:) ! Number of points on each latitude
logical :: luserpnm = .false. ! Use Belusov algorithm to compute RPNM array instead of per m
logical :: luseflt = .false. ! Use fast legendre transforms

! Extra inv_trans options
logical :: lvordiv = .false. ! Compute vorticity and divergence in grid point space
logical :: lscders = .false. ! Compute derivatives of scalar (North-South and East-West) in grid
                             ! point space
logical :: luvder = .false. ! Compute East-West derivatives of U and V wind in grid point space

! GSTATS options
logical :: lstats = .true. ! gstats statistics
logical :: ltrace_stats = .false.
logical :: lstats_omp = .false.
logical :: lstats_comms = .false.
logical :: lbarrier_stats = .false.
logical :: lbarrier_stats2 = .false.
logical :: ldetailed_stats = .false.
logical :: lstats_alloc = .false.
logical :: lsyncstats = .false.
logical :: lstatscpu = .false.
logical :: lstats_mem = .false.
logical :: lxml_stats = .false.
integer(kind=jpim) :: nstats_mem = 0
integer(kind=jpim) :: ntrace_stats = 0
integer(kind=jpim) :: nprnt_stats = 1
integer(kind=jpim) :: nopt_mem_tr = 0

logical :: lprint_norms = .false. ! Calculate and print spectral norms
logical :: lmeminfo = .false. ! Show information from FIAT routine ec_meminfo at the end

! The multiplier of the machine epsilon used as a tolerance for correctness checking
! ncheck = 0 (the default) means that correctness checking is disabled
integer(kind=jpim) :: ncheck = 0

logical :: lmpoff = .false. ! Message passing switch

! Verbosity level (0 or 1)
integer :: verbosity = 0

integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nthread
integer(kind=jpim) :: nprgpns ! Grid-point decomp
integer(kind=jpim) :: nprgpew ! Grid-point decomp
integer(kind=jpim) :: nprtrv = 0 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
integer(kind=jpim) :: mysetv
integer(kind=jpim) :: mysetw
integer(kind=jpim) :: mp_type = 2 ! Message passing type
integer(kind=jpim) :: mbx_size = 150000000 ! Mailbox size

integer(kind=jpim), allocatable :: numll(:), ivset(:), ivsetsc(:)
integer(kind=jpim) :: ivsetsc2(1)

integer(kind=jpim) :: nflevl

! sumpini
integer(kind=jpim) :: isqr
logical :: lsync_trans = .true. ! Activate barrier sync
logical :: leq_regions = .true. ! Eq regions flag

integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: ngpblks
! locals
integer(kind=jpim) :: iprtrv
integer(kind=jpim) :: iprtrw
integer(kind=jpim) :: iprused, ilevpp, irest, ilev, jlev

logical :: ldump_values = .false.

integer, external :: ec_mpirank
logical :: luse_mpi = .true.

character(len=16) :: cgrid = ''

integer(kind=jpim) :: ierr
integer :: icall_mode = 1
integer :: inum_wind_fields, inum_sc_3d_fields, inum_sc_2d_fields, itotal_fields
integer :: ipgp_start, ipgp_end, ipgpuv_start, ipgpuv_end

real(kind=jprb), allocatable :: global_field(:,:)

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "gath_grid.h"
#include "specnorm.h"
#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"
#include "ec_meminfo.intfb.h"

!===================================================================================================

luse_mpi = detect_mpirun()

! Setup
call get_command_line_arguments(nsmax, cgrid, iters, iters_warmup, nfld, nlev, lvordiv, lscders, luvder, &
  & luseflt, nopt_mem_tr, nproma, verbosity, ldump_values, lprint_norms, lmeminfo, nprtrv, nprtrw, ncheck, &
  & icall_mode)
if (cgrid == '') cgrid = cubic_octahedral_gaussian_grid(nsmax)
call parse_grid(cgrid, ndgl, nloen)
nflevg = nlev

!===================================================================================================

if (luse_mpi) then
  call mpl_init(ldinfo=(verbosity>=1))
  nproc  = mpl_nproc()
  myproc = mpl_myrank()
else
  nproc = 1
  myproc = 1
  mpl_comm = -1
  lsync_trans = .false.
endif
nthread = oml_max_threads()

call dr_hook_init()

!===================================================================================================

if( lstats ) call gstats(0,0)
ztinit = timef()

! only output to stdout on pe 1
if (nproc > 1) then
  if (myproc /= 1) then
    open(unit=nout, file='/dev/null')
  endif
endif

if (ldetailed_stats) then
  lstats_omp    = .true.
  lstats_comms  = .true.
  lstatscpu     = .true.
  nprnt_stats   = nproc
endif

!===================================================================================================

allocate(nprcids(nproc))
do jj = 1, nproc
  nprcids(jj) = jj
enddo

if (nproc <= 1) then
  lmpoff = .true.
endif

! Compute nprgpns and nprgpew
! This version selects most square-like distribution
! These will change if leq_regions=.true.
if (nproc == 0) nproc = 1
isqr = int(sqrt(real(nproc,jprb)))
do ja = isqr, nproc
  ib = nproc/ja
  if (ja*ib == nproc) then
    nprgpns = max(ja,ib)
    nprgpew = min(ja,ib)
    exit
  endif
enddo

! Compute nprtrv and nprtrw if not provided on the command line
if (nprtrv > 0 .or. nprtrw > 0) then
  if (nprtrv == 0) nprtrv = nproc/nprtrw
  if (nprtrw == 0) nprtrw = nproc/nprtrv
  if (nprtrw*nprtrv /= nproc) call abor1('ectrans_benchmark:nprtrw*nprtrv /= nproc')
else
  do jprtrv = 4, nproc
    nprtrv = jprtrv
    nprtrw = nproc/nprtrv
    if (nprtrv*nprtrw /= nproc) cycle
    if (nprtrv > nprtrw) exit
  enddo
  ! Go for approx square partition for backup
  if (nprtrv*nprtrw /= nproc .or. nprtrv > nprtrw) then
    isqr = int(sqrt(real(nproc,jprb)))
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

! Create communicators for mpi groups
if (.not.lmpoff) then
  call mpl_groups_create(nprtrw, nprtrv)
endif

if (lmpoff) then
  mysetw = (myproc - 1)/nprtrv + 1
  mysetv = mod(myproc - 1, nprtrv) + 1
else
  call mpl_cart_coords(myproc, mysetw, mysetv)

  ! Just checking for now...
  iprtrv = mod(myproc - 1, nprtrv) + 1
  iprtrw = (myproc - 1)/nprtrv + 1
  if (iprtrv /= mysetv .or. iprtrw /= mysetw) then
    call abor1('ectrans_benchmark:inconsistency when computing mysetw and mysetv')
  endif
endif

if (.not. lmpoff) then
  call mpl_buffer_method(kmp_type=mp_type, kmbx_size=mbx_size, kprocids=nprcids, ldinfo=(verbosity>=1))
endif

! Determine the number of levels attributed to each member of the V set
allocate(numll(nprtrv))
iprused = min(nflevg+1, nprtrv)
ilevpp = nflevg/nprtrv
irest = nflevg -ilevpp*nprtrv
do jroc = 1, nprtrv
  if (jroc <= irest) then
    numll(jroc) = ilevpp+1
  else
    numll(jroc) = ilevpp
  endif
enddo

nflevl = numll(mysetv)

ivsetsc2(1) = iprused
ifld = 0

!===================================================================================================
! Setup gstats
!===================================================================================================

if (lstats) then
  call gstats_setup(nproc, myproc, nprcids,                                            &
    & lstats, lstatscpu, lsyncstats, ldetailed_stats, lbarrier_stats, lbarrier_stats2, &
    & lstats_omp, lstats_comms, lstats_mem, nstats_mem, lstats_alloc,                  &
    & ltrace_stats, ntrace_stats, nprnt_stats, lxml_stats)
  call gstats_psut

  ! Assign labels to GSTATS regions
  call gstats_labels
endif

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

if (verbosity >= 1) write(nout,'(a)')'======= Setup ecTrans ======='

call gstats(1, 0)
call setup_trans0(kout=nout, kerr=nerr, kprintlev=merge(2, 0, verbosity == 1),                &
  &               kprgpns=nprgpns, kprgpew=nprgpew, kprtrw=nprtrw, ldsync_trans=lsync_trans,  &
  &               ldeq_regions=leq_regions, ldalloperm=.true., ldmpoff=.not.luse_mpi,         &
  &               kopt_memory_tr=nopt_mem_tr)
call gstats(1, 1)

call gstats(2, 0)
! IFS spectral fields are dimensioned NFLEVL, Nils !!
call set_ectrans_gpu_nflev(nflevl)
  ! We pass nflevl via environment variable in order not to change API
  ! In long run, ectrans should grow its internal buffers automatically
call setup_trans(ksmax=nsmax, kdgl=ndgl, kloen=nloen, ldsplit=.true., lduserpnm=luserpnm, &
  &              lduseflt=luseflt)
call gstats(2, 1)

call trans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg)

if (nproma == 0) then ! no blocking (default when not specified)
  nproma = ngptot
endif

! Calculate number of NPROMA blocks
ngpblks = (ngptot - 1)/nproma+1

!===================================================================================================
! Print information before starting
!===================================================================================================

! Print configuration details
if (verbosity >= 0 .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a)')'======= Start of runtime parameters ======='
  write(nout,'(" ")')
  write(nout,'("nsmax      ",i0)') nsmax
  write(nout,'("grid       ",a)') trim(cgrid)
  write(nout,'("ndgl       ",i0)') ndgl
  write(nout,'("nproc      ",i0)') nproc
  write(nout,'("nthread    ",i0)') nthread
  write(nout,'("nprgpns    ",i0)') nprgpns
  write(nout,'("nprgpew    ",i0)') nprgpew
  write(nout,'("nprtrw     ",i0)') nprtrw
  write(nout,'("nprtrv     ",i0)') nprtrv
  write(nout,'("ngptot     ",i0)') ngptot
  write(nout,'("ngptotg    ",i0)') ngptotg
  write(nout,'("nfld       ",i0)') nfld
  write(nout,'("nlev       ",i0)') nlev
  write(nout,'("nproma     ",i0)') nproma
  write(nout,'("ngpblks    ",i0)') ngpblks
  write(nout,'("nspec2     ",i0)') nspec2
  write(nout,'("nspec2g    ",i0)') nspec2g
  write(nout,'("luseflt    ",l1)') luseflt
  write(nout,'("nopt_mem_tr",i0)') nopt_mem_tr
  write(nout,'("lvordiv    ",l1)') lvordiv
  write(nout,'("lscders    ",l1)') lscders
  write(nout,'("luvder     ",l1)') luvder
  write(nout,'(" ")')
  write(nout,'(a)') '======= End of runtime parameters ======='
  write(nout,'(" ")')
end if

!===================================================================================================
! Allocate and initialize spectral arrays
!===================================================================================================

! Initialize vorticity and divergence - same for both call modes
allocate(zspvor(nflevl,nspec2))
allocate(zspdiv(nflevl,nspec2))
call initialize_spectral_field(nsmax, zspvor)
call initialize_spectral_field(nsmax, zspdiv)

! Initialize spectral arrays differently depending on call mode
if (icall_mode == 1) then
  allocate(zspscalar(nfld*nflevl+1,nspec2))
  call initialize_spectral_field(nsmax, zspscalar)
else
  allocate(zspsc3a(nflevl,nspec2,nfld))
  allocate(zspsc2(1,nspec2))
  do i = 1, nfld
    call initialize_spectral_field(nsmax, zspsc3a(:,:,i))
  enddo
  call initialize_spectral_field(nsmax, zspsc2)
endif

! Compute spectral distribution variables
allocate(ivset(nflevg))
ilev = 0
do jb = 1, nprtrv
  do jlev=1, numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo

allocate(ivsetsc(nfld*nflevg+1))
do i = 1, nfld
  ilev = 0
  do jb = 1, nprtrv
    do jlev = 1, numll(jb)
      ilev = ilev + 1
      ivsetsc(ilev + (i - 1)*nflevg) = jb
    enddo
  enddo
enddo
ivsetsc(nfld*nflevg+1) = 1

!===================================================================================================
! Allocate gridpoint arrays
!===================================================================================================

! Determine start and end slice points for grid point arrays when they are passed back to dir_trans
ipgp_start = 1
ipgp_end = (2 + nfld) * nflevg + 1
ipgpuv_start = 1
ipgpuv_end = 2

! Also enable vorticity divergence?
if (lvordiv) then
  inum_wind_fields = 4 ! Four fields - U, V, vorticity, divergence
  ! If lvordiv, skip the vorticity and divergence elements when passing zgp
  ! These two come first when enabled
  ipgp_start = ipgp_start + 2 * nflevg
  ipgp_end = ipgp_end + 2 * nflevg
  ipgpuv_start = ipgpuv_start + 2
  ipgpuv_end = ipgpuv_end + 2
else
  ! Otherwise just U and V
  inum_wind_fields = 2
endif

! Also make room for East-West derivatives of winds?
if (luvder) inum_wind_fields = inum_wind_fields + 2

! We always have our nfld 3D scalar fields
inum_sc_3d_fields = nfld

! We always have one 2D scalar field
inum_sc_2d_fields = 1

! Also make room for North-South and East-West derivatives of scalar fields
if (lscders) then
  inum_sc_3d_fields = inum_sc_3d_fields * 3
  inum_sc_2d_fields = inum_sc_2d_fields * 3
endif

! Finally, allocate grid point arrays
if (icall_mode == 1) then
  itotal_fields = nflevg * (inum_wind_fields + inum_sc_3d_fields) + inum_sc_2d_fields
  allocate(zgp(nproma,itotal_fields,ngpblks))
else
  allocate(zgpuv(nproma,nflevg,inum_wind_fields,ngpblks))
  allocate(zgp3a(nproma,nflevg,inum_sc_3d_fields,ngpblks))
  allocate(zgp2(nproma,inum_sc_2d_fields,ngpblks))
endif

!===================================================================================================
! Allocate norm arrays
!===================================================================================================

if (lprint_norms .or. ncheck > 0) then
  allocate(znormvor(nflevg))
  allocate(znormvor1(nflevg))
  allocate(znormdiv(nflevg))
  allocate(znormdiv1(nflevg))

  call specnorm(pspec=zspvor(1:nflevl,:), pnorm=znormvor1, kvset=ivset)
  call specnorm(pspec=zspdiv(1:nflevl,:), pnorm=znormdiv1, kvset=ivset)

  if (icall_mode == 1) then
    allocate(znormscalar(nfld*nflevg+1))
    allocate(znormscalar1(nfld*nflevg+1))
    call specnorm(pspec=zspscalar(:,:), pnorm=znormscalar1, kvset=ivsetsc)
  else
    allocate(znormsc3a(nflevg))
    allocate(znormsc3a1(nflevg))
    allocate(znormsc2(1))
    allocate(znormsc21(1))
    if (nfld > 0) call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormsc3a1, kvset=ivset)
    call specnorm(pspec=zspsc2(1:1,:), pnorm=znormsc21, kvset=ivsetsc2)
  endif

  if (verbosity >= 1 .and. myproc == 1) then
    do ifld = 1, nflevg
      write(nout,'("norm zspvor( ",i4,",:)   = ",f20.15)') ifld, znormvor1(ifld)
      write(nout,'("0x",Z16.16)') znormvor1(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv1(ifld)
      write(nout,'("0x",Z16.16)') znormdiv1(ifld)
    enddo
    if (icall_mode == 1) then
      do ifld = 1, nfld*nflevg+1
        write(nout,'("norm zspscalar(",i4,",:,1) = ",f20.15)') ifld, znormscalar1(ifld)
        write(nout,'("0x",Z16.16)') znormscalar1(ifld)
      enddo
    else
      if (nfld > 0) then
        do ifld = 1, nflevg
          write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormsc3a1(ifld)
          write(nout,'("0x",Z16.16)') znormsc3a1(ifld)
        enddo
      endif
      write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15)') 1, znormsc21(1)
      write(nout,'("0x",Z16.16)') znormsc21(1)
    endif
  endif
endif

!===================================================================================================
! Setup timers
!===================================================================================================

ztinit = (timef() - ztinit)/1000.0_jprd

if (verbosity >= 0 .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a,i0,a,f9.2,a)') "ectrans_benchmark initialisation, on ",nproc,&
                                & " tasks, took ",ztinit," sec"
  write(nout,'(" ")')
endif

if (iters <= 0) call abor1('ectrans_benchmark:iters <= 0')

allocate(ztstep(iters+iters_warmup))
allocate(ztstep1(iters+iters_warmup))
allocate(ztstep2(iters+iters_warmup))

if (verbosity >= 1 .and. myproc == 1) then
  write(nout,'(a)') '======= Start of spectral transforms  ======='
  write(nout,'(" ")')
endif



!===================================================================================================
! Do spectral transform loop
!===================================================================================================

gstats_lstats = .false.

write(nout,'(a,i0,a,i0,a)') 'Running for ', iters, ' iterations with ', iters_warmup, &
  & ' extra warm-up iterations'
write(nout,'(" ")')

do jstep = 1, iters+iters_warmup
  if (jstep == iters_warmup + 1) then
    gstats_lstats = .true.
    ztloop = timef()
  endif

  call gstats(3,0)
  ztstep(jstep) = timef()

  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

  ztstep1(jstep) = timef()
  call gstats(4,0)
  if (icall_mode == 1) then
    call inv_trans(pspvor=zspvor, pspdiv=zspdiv, pspscalar=zspscalar, pgp=zgp, &
      &            kvsetuv=ivset, kvsetsc=ivsetsc, &
      &            ldscders=lscders, ldvorgp=lvordiv, lddivgp=lvordiv, lduvder=luvder, &
      &            kproma=nproma)
  else
    call inv_trans(pspvor=zspvor, pspdiv=zspdiv, pspsc3a=zspsc3a, pspsc2=zspsc2, pgpuv=zgpuv, &
      &            pgp3a=zgp3a, pgp2=zgp2, &
      &            kvsetuv=ivset, kvsetsc2=ivsetsc2, kvsetsc3a=ivset, &
      &            ldscders=lscders, ldvorgp=lvordiv, lddivgp=lvordiv, lduvder=luvder, kproma=nproma)
  endif
  call gstats(4,1)

  ztstep1(jstep) = (timef() - ztstep1(jstep))/1000.0_jprd

  !=================================================================================================
  ! While in grid point space, dump the values to disk, for debugging only
  !=================================================================================================

  if (ldump_values .and. mod(jstep,10) == 1) then
    ! dump a field to a binary file
    if (myproc == 1) then
      allocate(global_field(ngptotg,1))
    endif
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgpuv(:,nflevg:nflevg,1,:), 'U', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgpuv(:,nflevg:nflevg,2,:), 'V', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp2(:,1:1,:), 'S', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp3a(:,nflevg:nflevg,1,:), 'T', noutdump)
    if (myproc == 1) then
      deallocate(global_field)
    endif
  endif

  !=================================================================================================
  ! Do direct transform
  !=================================================================================================

  ztstep2(jstep) = timef()

  call gstats(5,0)
  if (icall_mode == 1) then
    call dir_trans(pgp=zgp(:,ipgp_start:ipgp_end,:), pspvor=zspvor, pspdiv=zspdiv, &
      &            pspscalar=zspscalar, kvsetuv=ivset, kvsetsc=ivsetsc, kproma=nproma)
  else
    call dir_trans(pgpuv=zgpuv(:,:,ipgpuv_start:ipgpuv_end,:), &
      &            pgp3a=zgp3a(:,:,1:nfld,:), pgp2=zgp2(:,1:1,:), &
      &            pspvor=zspvor, pspdiv=zspdiv, pspsc3a=zspsc3a, pspsc2=zspsc2, &
      &            kvsetuv=ivset, kvsetsc2=ivsetsc2, kvsetsc3a=ivset, kproma=nproma)
  endif
  call gstats(5,1)
  ztstep2(jstep) = (timef() - ztstep2(jstep))/1000.0_jprd

  ztstep(jstep) = (timef() - ztstep(jstep))/1000.0_jprd

  !=================================================================================================
  ! Print norms
  !=================================================================================================

  if (lprint_norms) then
    call gstats(6,0)
    call specnorm(pspec=zspvor(1:nflevl,:), pnorm=znormvor, kvset=ivset)
    call specnorm(pspec=zspdiv(1:nflevl,:), pnorm=znormdiv, kvset=ivset)

    if (icall_mode == 1) then
      call specnorm(pspec=zspscalar(:,:), pnorm=znormscalar, kvset=ivsetsc)
    else
      if (nfld > 0) call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormsc3a, kvset=ivset)
      call specnorm(pspec=zspsc2(1:1,:), pnorm=znormsc2, kvset=ivsetsc2)
    endif

    if (myproc == 1) then
      zmaxerr(1) = maxval(abs((znormvor1 / znormvor) - 1.0_jprb))
      zmaxerr(2) = maxval(abs((znormdiv1 / znormdiv) - 1.0_jprb))
      if (icall_mode == 1) then
        zmaxerr(3) = maxval(abs((znormscalar1 / znormscalar) - 1.0_jprb))
        write(nout,'("time step ",i6," took", f8.4," | zspvor max err="e10.3,&
        & " | zspdiv max err="e10.3," | zspscalar max err="e10.3)') &
        &  jstep, ztstep(jstep), zmaxerr(1), zmaxerr(2), zmaxerr(3)
      else
        zmaxerr(4) = maxval(abs((znormsc21 / znormsc2) - 1.0_jprb))
        if (nfld > 0) then
          zmaxerr(3) = maxval(abs((znormsc3a1 / znormsc3a) - 1.0_jprb))
          write(nout,'("time step ",i6," took", f8.4," | zspvor max err="e10.3,&
          & " | zspdiv max err="e10.3," | zspsc3a max err="e10.3," | zspsc2 max err="e10.3)') &
          &  jstep, ztstep(jstep), zmaxerr(1), zmaxerr(2), zmaxerr(3), zmaxerr(4)
        else
          write(nout,'("time step ",i6," took", f8.4," | zspvor max err="e10.3,&
                      & " | zspdiv max err="e10.3," | zspsc2 max err="e10.3)') &
                      &  jstep, ztstep(jstep), zmaxerr(1), zmaxerr(2), zmaxerr(4)
        endif
      endif
    endif
    call gstats(6,1)
  else
    write(nout,'("Time step ",i6," took", f8.4)') jstep, ztstep(jstep)
  endif
  call gstats(3,1)
enddo

!===================================================================================================

ztloop = (timef() - ztloop)/1000.0_jprd

write(nout,'(" ")')
write(nout,'(a)') '======= End of spectral transforms  ======='
write(nout,'(" ")')

if (lprint_norms .or. ncheck > 0) then
  call specnorm(pspec=zspvor(1:nflevl,:), pnorm=znormvor, kvset=ivset)
  call specnorm(pspec=zspdiv(1:nflevl,:), pnorm=znormdiv, kvset=ivset)

  if (icall_mode == 1) then
    call specnorm(pspec=zspscalar(:,:), pnorm=znormscalar, kvset=ivsetsc)
  else
    if (nfld > 0) call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormsc3a, kvset=ivset)
    call specnorm(pspec=zspsc2(1:1,:), pnorm=znormsc2, kvset=ivsetsc2)
  endif

  if (myproc == 1) then
    zmaxerr = -99.0_jprd
    zmaxerr(1) = maxval(abs((real(znormvor1,jprd) / (real(znormvor,jprd)) - 1.0_jprd)))
    if (verbosity >= 1) then
      do ifld = 1, nflevg
        write(nout,'("norm zspvor( ",i4,")     = ",f20.15)') ifld, znormvor(ifld)
        write(nout,'("0x",Z16.16)') znormvor(ifld)
      enddo
    endif
    zmaxerr(2) = maxval(abs((real(znormdiv1,jprd) / (real(znormdiv,jprd)) - 1.0_jprd)))
    if (verbosity >= 1) then
      do ifld = 1, nflevg
        write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv(ifld)
        write(nout,'("0x",Z16.16)') znormdiv(ifld)
      enddo
    endif
    if (icall_mode == 1) then
      zmaxerr(3) = maxval(abs((znormscalar1 / znormscalar) - 1.0_jprb))
      if (verbosity >= 1) then
        do ifld = 1, nfld*nflevg+1
          write(nout,'("norm znormscalar( ",i4,",:)   = ",f20.15)') ifld, znormscalar(ifld)
          write(nout,'("0x",Z16.16)') znormscalar(ifld)
        enddo
      endif
    else
      zmaxerr(4) = maxval(abs((znormsc21 / znormsc2) - 1.0_jprb))
      if (verbosity >= 1) then
        write(nout,'("norm znormsc2( ",i4,",:)   = ",f20.15)') 1, znormsc2(1)
        write(nout,'("0x",Z16.16)') znormsc2(1)
      endif
      if (nfld > 0) then
        zmaxerr(3) = maxval(abs((znormsc3a1 / znormsc3a) - 1.0_jprb))
        if (verbosity >= 1) then
          do ifld = 1, nflevg
            write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormsc3a(ifld)
            write(nout,'("0x",Z16.16)') znormsc3a(ifld)
          enddo
        endif
      endif
    endif

    ! maximum error across all fields
    zmaxerrg = maxval(zmaxerr)

    if (verbosity >= 1) write(nout,*)
    write(nout,'("max error zspvor(1:nlev,:)    = ",e10.3)') zmaxerr(1)
    write(nout,'("max error zspdiv(1:nlev,:)    = ",e10.3)') zmaxerr(2)
    if (icall_mode == 1) then
      write(nout,'("max error zspscalar(1:nlev,:,1) = ",e10.3)') zmaxerr(3)
    else
      if (nfld > 0) write(nout,'("max error zspsc3a(1:nlev,:,1) = ",e10.3)') zmaxerr(3)
      write(nout,'("max error zspsc2(1:1,:)       = ",e10.3)') zmaxerr(4)
    endif
    write(nout,*)
    write(nout,'("max error combined =          = ",e10.3)') zmaxerrg
    write(nout,*)
  endif
  if (ncheck > 0) then
    ierr = 0
    if (myproc == 1) then
      ! If the maximum spectral norm error across all fields is greater than 100 times the machine
      ! epsilon, fail the test
      if (zmaxerrg > real(ncheck, jprb) * epsilon(1.0_jprb)) then
        write(nout, '(a)') '*******************************'
        write(nout, '(a)') 'Correctness test failed'
        write(nout, '(a,1e9.2)') 'Maximum spectral norm error = ', zmaxerrg
        write(nout, '(a,1e9.2)') 'Error tolerance = ', real(ncheck, jprb) * epsilon(1.0_jprb)
        write(nout, '(a)') '*******************************'
        ierr = 1
      endif
    endif

    ! Root rank broadcasts the correctness checker result to the other ranks
    if (luse_mpi) then
      call mpl_broadcast(ierr,kroot=1,ktag=1)
    endif

    ! Halt if correctness checker failed
    if (ierr == 1) then
      error stop
    endif
  endif
endif

!===================================================================================================
! Calculate timings
!===================================================================================================

ztstepavg = sum(ztstep(iters_warmup+1:))
ztstepmin = minval(ztstep(iters_warmup+1:))
ztstepmax = maxval(ztstep(iters_warmup+1:))
ztstepavg1 = sum(ztstep1(iters_warmup+1:))
ztstepmin1 = minval(ztstep1(iters_warmup+1:))
ztstepmax1 = maxval(ztstep1(iters_warmup+1:))
ztstepavg2 = sum(ztstep2(iters_warmup+1:))
ztstepmin2 = minval(ztstep2(iters_warmup+1:))
ztstepmax2 = maxval(ztstep2(iters_warmup+1:))

if (luse_mpi) then
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
endif

ztstepavg = (ztstepavg/real(nproc,jprb))/real(iters,jprd)
ztloop = ztloop/real(nproc,jprd)
ztstep(:) = ztstep(:)/real(nproc,jprd)
ztstepmed = get_median(ztstep(iters_warmup+1:))

ztstepavg1 = (ztstepavg1/real(nproc,jprb))/real(iters,jprd)
ztstep1(:) = ztstep1(:)/real(nproc,jprd)
ztstepmed1 = get_median(ztstep1(iters_warmup+1:))

ztstepavg2 = (ztstepavg2/real(nproc,jprb))/real(iters,jprd)
ztstep2(:) = ztstep2(:)/real(nproc,jprd)
ztstepmed2 = get_median(ztstep2(iters_warmup+1:))

write(nout,'(a)') '======= Start of time step stats ======='
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
write(nout,'(a)') '======= End of time step stats ======='
write(nout,'(" ")')

if (lstack) then
  ! Gather stack usage statistics
  istack = getstackusage()
  if (myproc == 1) then
    print 9000, istack
    9000 format("Stack utilisation information",/,&
         &"=============================",//,&
         &"Task           size(bytes)",/,&
         &"====           ===========",//,&
         &"   1",11x,i10)

    do i = 2, nproc
      call mpl_recv(istack, ksource=nprcids(i), ktag=i, cdstring='ectrans_benchmark:')
      print '(i4,11x,i10)', i, istack
    enddo
  else
    call mpl_send(istack, kdest=nprcids(1), ktag=myproc, cdstring='ectrans_benchmark:')
  endif
endif


!===================================================================================================
! Cleanup
!===================================================================================================



!===================================================================================================

if (lstats) then
  call gstats(0,1)
  call gstats_print(nout, zaveave, jpmaxstat)
endif

if (lmeminfo) then
  write(nout,*)
  call ec_meminfo(nout, "", mpl_comm, kbarr=1, kiotask=-1, &
      & kcall=1)
endif

!===================================================================================================
! Finalize MPI
!===================================================================================================

if (luse_mpi) then
  call mpl_end(ldmeminfo=.false.)
endif

!===================================================================================================
! Close file
!===================================================================================================

if (nproc > 1) then
  if (myproc /= 1) then
    close(unit=nout)
  endif
endif

!===================================================================================================

contains

!===================================================================================================

subroutine parse_grid(cgrid,ndgl,nloen)

  character(len=*) :: cgrid
  integer, intent(inout) :: ndgl
  integer, intent(inout), allocatable :: nloen(:)
  integer :: ios
  integer :: gaussian_number
  read(cgrid(2:len_trim(cgrid)),*,IOSTAT=ios) gaussian_number
  if (ios==0) then
    ndgl = 2 * gaussian_number
    allocate(nloen(ndgl))
    if (cgrid(1:1) == 'F') then ! Regular Gaussian grid
      nloen(:) = gaussian_number * 4
      return
    endif
    if (cgrid(1:1) == 'O') then ! Octahedral Gaussian grid
      do i = 1, ndgl / 2
        nloen(i) = 20 + 4 * (i - 1)
        nloen(ndgl - i + 1) = nloen(i)
      end do
      return
    endif
  endif
  call parsing_failed("ERROR: Unsupported grid specified: "// trim(cgrid))

end subroutine

!===================================================================================================

subroutine str2int(str, int, stat)

  character(len=*), intent(in) :: str
  integer, intent(out) :: int
  integer, intent(out) :: stat
  read(str, *, iostat=stat) int

end subroutine str2int

!===================================================================================================

function get_int_value(cname, iarg) result(value)

  integer :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg
  character(len=128) :: carg
  integer :: stat

  carg = get_str_value(cname, iarg)
  call str2int(carg, value, stat)

  if (stat /= 0) then
    call parsing_failed("Invalid argument for " // trim(cname) // ": " // trim(carg))
  end if

end function

!===================================================================================================

function get_str_value(cname, iarg) result(value)

  character(len=128) :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg

  iarg = iarg + 1
  call get_command_argument(iarg, value)

  if (value == "") then
    call parsing_failed("Invalid argument for " // trim(cname) // ": no value provided")
  end if

end function

!===================================================================================================

subroutine print_help(unit)

  integer, optional :: unit
  integer :: nout = 6
  if (present(unit)) then
    nout = unit
  endif

  write(nout, "(a)") ""

  if (jprb == jprd) then
    write(nout, "(a)") "NAME    ectrans-benchmark-" // VERSION // "-dp"
  else
    write(nout, "(a)") "NAME    ectrans-benchmark-" // VERSION // "-sp"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "DESCRIPTION"
  write(nout, "(a)") "        This program tests ecTrans by transforming fields back and forth&
    & between spectral "
  if (jprb == jprd) then
    write(nout, "(a)") "        space and grid-point space (double-precision version)"
  else
    write(nout, "(a)") "        space and grid-point space (single-precision version)"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "USAGE"
  if (jprb == jprd) then
    write(nout, "(a)") "        ectrans-benchmark-" // VERSION // "-dp [options]"
  else
    write(nout, "(a)") "        ectrans-benchmark-" // VERSION // "-sp [options]"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "OPTIONS"
  write(nout, "(a)") "    -h, --help          Print this message"
  write(nout, "(a)") "    -v                  Run with verbose output"
  write(nout, "(a)") "    -t, --truncation T  Run with this triangular spectral truncation&
    & (default = 79)"
  write(nout, "(a)") "    -g, --grid GRID     Run with this grid. Possible values: O<N>, F<N>"
  write(nout, "(a)") "                        If not specified, O<N> is used with N=truncation+1&
    & (cubic relation)"
  write(nout, "(a)") "    -n, --niter NITER   Run for this many inverse/direct transform&
    & iterations (default = 10)"
  write(nout, "(a)") "    --niter-warmup      Number of warm up iterations,&
    & for which timing statistics should be ignored (default = 3)"
  write(nout, "(a)") "    -f, --nfld NFLD     Number of scalar fields (default = 1)"
  write(nout, "(a)") "    -l, --nlev NLEV     Number of vertical levels (default = 1)"
  write(nout, "(a)") "    --vordiv            Also transform vorticity-divergence to wind"
  write(nout, "(a)") "    --scders            Compute scalar derivatives (default off)"
  write(nout, "(a)") "    --uvders            Compute uv East-West derivatives (default off). Only&
    & when also --vordiv is given"
  write(nout, "(a)") "    --flt               Run with fast Legendre transforms (default off)"
  write(nout, "(a)") "    --nproma NPROMA     Run with NPROMA (default no blocking: NPROMA=ngptot)"
  write(nout, "(a)") "    --norms             Calculate and print spectral norms of transformed&
    & fields"
  write(nout, "(a)") "                        The computation of spectral norms will skew overall&
    & timings"
  write(nout, "(a)") "    --meminfo           Show diagnostic information from FIAT's ec_meminfo&
    & subroutine on memory usage, thread-binding etc."
  write(nout, "(a)") "    --nprtrv            Size of V set in spectral decomposition"
  write(nout, "(a)") "    --nprtrw            Size of W set in spectral decomposition"
  write(nout, "(a)") "    -c, --check VALUE   The multiplier of the machine epsilon used as a&
   & tolerance for correctness checking"
  write(nout, "(a)") "    --callmode          The call mode for INV_TRANS and DIR_TRANS (1 or 2)"
  write(nout, "(a)") ""
  write(nout, "(a)") "DEBUGGING"
  write(nout, "(a)") "    --dump-values       Output gridpoint fields in unformatted binary file"
  write(nout, "(a)") ""

end subroutine print_help

!===================================================================================================

subroutine parsing_failed(message)

  character(len=*), intent(in) :: message
  if (luse_mpi) call mpl_init(ldinfo=.false.)
  if (ec_mpirank() == 0) then
    write(nerr,"(a)") trim(message)
    call print_help(unit=nerr)
  endif
  if (luse_mpi) call mpl_end(ldmeminfo=.false.)
  stop

end subroutine

!===================================================================================================

subroutine get_command_line_arguments(nsmax, cgrid, iters, iters_warmup, nfld, nlev, lvordiv, lscders, luvder, &
  &                                   luseflt, nopt_mem_tr, nproma, verbosity, ldump_values, lprint_norms, &
  &                                   lmeminfo, nprtrv, nprtrw, ncheck, icall_mode)

#ifdef _OPENACC
  use openacc, only: acc_init, acc_get_device_type
#endif

  integer, intent(inout) :: nsmax           ! Spectral truncation
  character(len=16), intent(inout) :: cgrid ! Spectral truncation
  integer, intent(inout) :: iters           ! Number of iterations for transform test
  integer, intent(inout) :: iters_warmup    ! Number of iterations for transform test
  integer, intent(inout) :: nfld            ! Number of scalar fields
  integer, intent(inout) :: nlev            ! Number of vertical levels
  logical, intent(inout) :: lvordiv         ! Also transform vorticity/divergence
  logical, intent(inout) :: lscders         ! Compute scalar derivatives
  logical, intent(inout) :: luvder          ! Compute uv East-West derivatives
  logical, intent(inout) :: luseflt         ! Use fast Legendre transforms
  integer, intent(inout) :: nopt_mem_tr     ! Use of heap or stack memory for ZCOMBUF arrays in transposition arrays (0 for heap, 1 for stack)
  integer, intent(inout) :: nproma          ! NPROMA
  integer, intent(inout) :: verbosity       ! Level of verbosity
  logical, intent(inout) :: ldump_values    ! Dump values of grid point fields for debugging
  logical, intent(inout) :: lprint_norms    ! Calculate and print spectral norms of fields
  logical, intent(inout) :: lmeminfo        ! Show information from FIAT ec_meminfo routine at the
                                            ! end
  integer, intent(inout) :: nprtrv          ! Size of V set (spectral decomposition)
  integer, intent(inout) :: nprtrw          ! Size of W set (spectral decomposition)
  integer, intent(inout) :: ncheck          ! The multiplier of the machine epsilon used as a
                                            ! tolerance for correctness checking
  integer, intent(inout) :: icall_mode      ! The call mode for inv_trans and dir_trans
                                            ! 1: pspvor, pspdiv, pspscalar, pgp
                                            ! 2: pspvor, pspdiv, pspsc3a, pspsc2, pgpuv, pgp3a, pgp2

  character(len=128) :: carg          ! Storage variable for command line arguments
  integer            :: iarg = 1      ! Argument index

#ifdef _OPENACC
  call acc_init(acc_get_device_type())
#endif

  do while (iarg <= command_argument_count())
    call get_command_argument(iarg, carg)

    select case(carg)
      ! Parse help argument
      case('-h', '--help')
        if (luse_mpi) call mpl_init(ldinfo=.false.)
        if (ec_mpirank()==0) call print_help()
        if (luse_mpi) call mpl_end(ldmeminfo=.false.)
        stop
      ! Parse verbosity argument
      case('-v')
        verbosity = 1
      ! Parse number of iterations argument
      case('-n', '--niter')
        iters = get_int_value('-n', iarg)
        if (iters < 1) then
          call parsing_failed("Invalid argument for -n: must be > 0")
        end if
      case('--niter-warmup')
        iters_warmup = get_int_value('--niter-warmup', iarg)
        if (iters_warmup < 0) then
          call parsing_failed("Invalid argument for --niter-warmup: must be >= 0")
        end if
      ! Parse spectral truncation argument
      case('-t', '--truncation')
        nsmax = get_int_value('-t', iarg)
        if (nsmax < 1) then
          call parsing_failed("Invalid argument for -t: must be > 0")
        end if
      case('-g', '--grid'); cgrid = get_str_value('-g', iarg)
      case('-f', '--nfld'); nfld = get_int_value('-f', iarg)
      case('-l', '--nlev'); nlev = get_int_value('-l', iarg)
      case('--vordiv'); lvordiv = .True.
      case('--scders'); lscders = .True.
      case('--uvders'); luvder = .True.
      case('--flt'); luseflt = .True.
      case('--mem-tr'); nopt_mem_tr = get_int_value('--mem-tr', iarg)
      case('--nproma'); nproma = get_int_value('--nproma', iarg)
      case('--dump-values'); ldump_values = .true.
      case('--norms'); lprint_norms = .true.
      case('--meminfo'); lmeminfo = .true.
      case('--nprtrv'); nprtrv = get_int_value('--nprtrv', iarg)
      case('--nprtrw'); nprtrw = get_int_value('--nprtrw', iarg)
      case('-c', '--check'); ncheck = get_int_value('-c', iarg)
      case('--callmode')
          icall_mode = get_int_value('--callmode', iarg)
          if (icall_mode /= 1 .and. icall_mode /= 2) then
            call parsing_failed("Invalid argument for --callmode: must be 1 or 2")
          end if
      case default
        call parsing_failed("Unrecognised argument: " // trim(carg))

    end select
    iarg = iarg + 1
  end do

  if (.not. lvordiv) then
    luvder = .false.
  endif

end subroutine get_command_line_arguments
!===================================================================================================

function cubic_octahedral_gaussian_grid(nsmax) result(cgrid)

  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'O',nsmax+1

end function


!===================================================================================================

function get_median(vec) result(median)

  real(kind=jprd), intent(in) :: vec(:)
  real(kind=jprd) :: median

  real(kind=jprd) :: vec_sorted(size(vec))
  real(kind=jprd) :: x

  integer :: i, j, n

  n = size(vec)

  ! Sort in ascending order
  vec_sorted = vec
  do i = 2, n
    x = vec_sorted(i)
    j = i - 1
    do while (j >= 1)
      if (vec_sorted(j) <= x) exit
      vec_sorted(j + 1) = vec_sorted(j)
      j = j - 1
    end do
    vec_sorted(j + 1) = x
  end do

  ! Calculate median according to if there is an even or odd number of elements
  if (mod(n, 2) == 0) then
    median = (vec_sorted(n/2) + vec_sorted(n/2+1))/2.0_jprd
  else
    median = vec_sorted((n+1)/2)
  endif

end function get_median

!===================================================================================================

subroutine initialize_spectral_field(nsmax, field)

  integer,         intent(in)    :: nsmax      ! Spectral truncation
  real(kind=jprb), intent(inout) :: field(:,:) ! Field to initialize

  integer :: i

  do i = 1, size(field,1)
    call initialize_2d_spectral_field(nsmax, field(i,:))
  enddo

end subroutine initialize_spectral_field

!===================================================================================================

subroutine initialize_2d_spectral_field(nsmax, field)

  integer,         intent(in)    :: nsmax    ! Spectral truncation
  real(kind=jprb), intent(inout) :: field(:) ! Field to initialize

  integer :: index, num_my_zon_wns
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

subroutine dump_gridpoint_field(jstep, myproc, nproma, gfld, fld, fldchar, noutdump)

  ! Dump a 2d field to a binary file.

  integer(kind=jpim), intent(in) :: jstep ! Time step, used for naming file
  integer(kind=jpim), intent(in) :: myproc ! MPI rank, used for naming file
  integer(kind=jpim), intent(in) :: nproma ! Size of nproma
  real(kind=jprb)   , intent(inout) :: gfld(:,:) ! 2d global field
  real(kind=jprb)   , intent(in) :: fld(:,:,:) ! 3d local field
  character         , intent(in) :: fldchar ! Single character field identifier
  integer(kind=jpim), intent(in) :: noutdump ! Tnit number for output file

  character(len=10) :: filename = "x.xxxx.dat"

  if (myproc == 1) then
    write(filename(1:1),'(a1)') fldchar
    write(filename(3:6),'(i4.4)') jstep
    open(noutdump,file=filename,form='unformatted')
  endif
  do ilev=1,size(fld,2)
    call gath_grid(gfld(:,:),nproma,1,(/1/),1,fld(:,ilev:ilev,:))
    if (myproc == 1) write(unit=noutdump) gfld(:,1)
  enddo
  if (myproc == 1) then
    close(noutdump)
  endif

end subroutine dump_gridpoint_field

!===================================================================================================

function detect_mpirun() result(lmpi_required)
  use ec_env_mod, only : ec_putenv
  logical :: lmpi_required
  integer :: ilen
  integer, parameter :: nvars = 4
  character(len=32), dimension(nvars) :: cmpirun_detect
  character(len=4) :: clenv
  integer :: ivar

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  cmpirun_detect(1) = 'OMPI_COMM_WORLD_SIZE'  ! openmpi
  cmpirun_detect(2) = 'ALPS_APP_PE'           ! cray pe
  cmpirun_detect(3) = 'PMI_SIZE'              ! intel
  cmpirun_detect(4) = 'SLURM_NTASKS'          ! slurm

  lmpi_required = .false.
  do ivar = 1, nvars
    call get_environment_variable(name=trim(cmpirun_detect(ivar)), length=ilen)
    if (ilen > 0) then
      lmpi_required = .true.
      exit ! break
    endif
  enddo

  call get_environment_variable(name="ECTRANS_USE_MPI", value=clenv, length=ilen )
  if (ilen > 0) then
      lmpi_required = .true.
      if( trim(clenv) == "0" .or. trim(clenv) == "OFF" .or. trim(CLENV) == "off" .or. trim(clenv) == "F" ) then
        lmpi_required = .false.
      endif
      call ec_putenv("DR_HOOK_ASSERT_MPI_INITIALIZED=0", overwrite=.true.)
  endif
end function

!===================================================================================================

! Assign GSTATS labels to the main regions of ecTrans
subroutine gstats_labels

  call gstats_label(0,   '   ', 'PROGRAM        - Total')
  call gstats_label(1,   '   ', 'SETUP_TRANS0   - Setup ecTrans')
  call gstats_label(2,   '   ', 'SETUP_TRANS    - Setup ecTrans handle')
  call gstats_label(3,   '   ', 'TIME STEP      - Time step')
  call gstats_label(4,   '   ', 'INV_TRANS      - Inverse transform')
  call gstats_label(5,   '   ', 'DIR_TRANS      - Direct transform')
  call gstats_label(6,   '   ', 'NORMS          - Norm comp. (optional)')
  call gstats_label(102, '   ', 'LTINV_CTL      - Inv. Legendre transform')
  call gstats_label(103, '   ', 'LTDIR_CTL      - Dir. Legendre transform')
  call gstats_label(106, '   ', 'FTDIR_CTL      - Dir. Fourier transform')
  call gstats_label(107, '   ', 'FTINV_CTL      - Inv. Fourier transform')
  call gstats_label(140, '   ', 'SULEG          - Comp. of Leg. poly.')
  call gstats_label(152, '   ', 'LTINV_CTL      - M to L transposition')
  call gstats_label(153, '   ', 'LTDIR_CTL      - L to M transposition')
  call gstats_label(157, '   ', 'FTINV_CTL      - L to G transposition')
  call gstats_label(158, '   ', 'FTDIR_CTL      - G to L transposition')
  call gstats_label(400, '   ', 'GSTATS         - GSTATS itself')

end subroutine gstats_labels

!===================================================================================================

subroutine set_ectrans_gpu_nflev(kflev)
  use ec_env_mod, only : ec_putenv
  integer(kind=jpim), intent(in) :: kflev
  character(len=32) :: ECTRANS_GPU_NFLEV
  write(ECTRANS_GPU_NFLEV,'(A,I0)') "ECTRANS_GPU_NFLEV=",kflev
  call ec_putenv(ECTRANS_GPU_NFLEV, overwrite=.true.)
end subroutine

end program ectrans_benchmark

!===================================================================================================
