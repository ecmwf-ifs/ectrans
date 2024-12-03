program ectrans_lam_benchmark

!
! Spectral transform test for Limited-Area geometry
!
! This test performs spectral to real and real to spectral transforms repeated in
! timed loop.
!
! 1) One "surface" field is always transformed:
!      zspsc2(1,1:nspec2) <-> zgmvs(1:nproma,1:1,1:ngbplk)
!
! 2) A Multiple "3d" fields are transformed and can be disabled with "--nfld 0"
!
!      zspsc3a(1:nlev,1:nspec2,1:nfld) <-> zgp3a(1:nproma,1:nlev,1:nfld,1:ngpblk)
!
! 3) Optionally a "3d" vorticity/divergence field is transformed to uv (wind) and
!   can be enabled with "--vordiv"
!
!      zspvor(1:nlev,1:nspec2) / zspdiv(1:nlev,1:nspec2) <-> zgpuv(1:nproma,1:nlev,1:2,1:ngpblk)
!
! 4) Optionally scalar derivatives can be computed for the fields described in 1) and 2)
!    This must be enabled with "--scders"
!
! 5) Optionally uv East-West derivate can be computed from vorticity/divergence.
!    This must be enabled with "--vordiv --uvders"
!
!
! Authors : George Mozdzynski
!           Willem Deconinck
!           Ioan Hadade
!           Sam Hatfield
!           Daan Degrauwe

use parkind1, only: jpim, jprb, jprd
use oml_mod ,only : oml_max_threads
use omp_lib, only: omp_get_wtime
use mpl_module
use yomgstats, only: jpmaxstat
use yomhook, only : dr_hook_init

implicit none

integer(kind=jpim) :: istack, getstackusage
real(kind=jprb), dimension(1) :: zmaxerr(5), zerr(5)
real(kind=jprb) :: zmaxerrg

! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT
integer(kind=jpim), parameter :: noutdump = 7 ! Unit number for field output

! Default parameters
integer(kind=jpim) :: nlon    = 128   ! Zonal dimension
integer(kind=jpim) :: nlat    = 128   ! Meridional dimension
integer(kind=jpim) :: nsmax   = 0   ! Spectral meridional truncation
integer(kind=jpim) :: nmsmax  = 0   ! Spectral zonal truncation
integer(kind=jpim) :: iters   = 10  ! Number of iterations for transform test
integer(kind=jpim) :: nfld    = 1   ! Number of scalar fields 
integer(kind=jpim) :: nlev    = 1   ! Number of vertical levels

integer(kind=jpim) :: nloen(1)  ! only one value needed for LAM
integer(kind=jpim) :: nflevg
integer(kind=jpim) :: nspec2
integer(kind=jpim) :: ngptot
integer(kind=jpim) :: ngptotg
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

real(kind=jprd) :: ztinit, ztloop, ztstepmax, ztstepmin, ztstepavg, ztstepmed
real(kind=jprd) :: ztstepmax1, ztstepmin1, ztstepavg1, ztstepmed1
real(kind=jprd) :: ztstepmax2, ztstepmin2, ztstepavg2, ztstepmed2
real(kind=jprd), allocatable :: ztstep(:), ztstep1(:), ztstep2(:)

real(kind=jprb), allocatable :: znormsp(:), znormsp0(:), znormdiv(:), znormdiv0(:)
real(kind=jprb), allocatable :: znormvor(:), znormvor0(:), znormt(:), znormt0(:)
real(kind=jprd) :: zaveave(0:jpmaxstat)

! Grid-point space data structures
real(kind=jprb), allocatable, target :: zgmv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), allocatable, target :: zgmvs  (:,:,:)   ! Single level fields at t and t-dt
real(kind=jprb), pointer :: zgp3a (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgpuv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgp2 (:,:,:) ! Single level fields at t and t-dt

! Spectral space data structures
real(kind=jprb), allocatable, target :: sp3d(:,:,:)
real(kind=jprb), pointer :: zspvor(:,:) => null()
real(kind=jprb), pointer :: zspdiv(:,:) => null()
real(kind=jprb), pointer :: zspsc3a(:,:,:) => null()
real(kind=jprb), allocatable :: zspsc2(:,:)
real(kind=jprb), allocatable :: zmeanu(:), zmeanv(:)

logical :: lstack = .false. ! Output stack info
logical :: luserpnm = .false.
logical :: lkeeprpnm = .false.
logical :: ltrace_stats = .false.
logical :: lstats_omp = .false.
logical :: lstats_comms = .false.
logical :: lstats_mpl = .false.
logical :: lstats = .false. ! gstats statistics
logical :: lbarrier_stats = .false.
logical :: lbarrier_stats2 = .false.
logical :: ldetailed_stats = .false.
logical :: lstats_alloc = .false.
logical :: lsyncstats = .false.
logical :: lstatscpu = .false.
logical :: lstats_mem = .false.
logical :: lxml_stats = .false.
logical :: lfftw = .false. ! Use FFTW for Fourier transforms
logical :: lvordiv = .false.
logical :: lscders = .false.
logical :: luvders = .false.
logical :: lprint_norms = .false. ! Calculate and print spectral norms
logical :: lmeminfo = .false. ! Show information from FIAT routine ec_meminfo at the end

integer(kind=jpim) :: nstats_mem = 0
integer(kind=jpim) :: ntrace_stats = 0
integer(kind=jpim) :: nprnt_stats = 1

! The multiplier of the machine epsilon used as a tolerance for correctness checking
! ncheck = 0 (the default) means that correctness checking is disabled
integer(kind=jpim) :: ncheck = 0

logical :: lmpoff = .false. ! Message passing switch

! Verbosity level (0 or 1)
integer :: verbosity = 0

integer(kind=jpim) :: nmax_resol = 37 ! Max number of resolutions
integer(kind=jpim) :: npromatr = 0 ! nproma for trans lib

integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nthread
integer(kind=jpim) :: nprgpns = 0 ! Grid-point decomp
integer(kind=jpim) :: nprgpew = 0 ! Grid-point decomp
integer(kind=jpim) :: nprtrv = 0 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
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
logical :: lsync_trans = .false. ! Activate barrier sync


integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: ngpblks
! locals
integer(kind=jpim) :: iprtrv
integer(kind=jpim) :: iprtrw
integer(kind=jpim) :: iprused, ilevpp, irest, ilev, jlev

integer(kind=jpim) :: ndimgmv  = 0 ! Third dim. of gmv "(nproma,nflevg,ndimgmv,ngpblks)"
integer(kind=jpim) :: ndimgmvs = 0 ! Second dim. gmvs "(nproma,ndimgmvs,ngpblks)"

integer(kind=jpim) :: jbegin_uv = 0
integer(kind=jpim) :: jend_uv   = 0
integer(kind=jpim) :: jbegin_sc = 0
integer(kind=jpim) :: jend_sc   = 0
integer(kind=jpim) :: jbegin_scder_NS = 0
integer(kind=jpim) :: jend_scder_NS = 0
integer(kind=jpim) :: jbegin_scder_EW = 0
integer(kind=jpim) :: jend_scder_EW = 0
integer(kind=jpim) :: jbegin_uder_EW = 0
integer(kind=jpim) :: jend_uder_EW = 0
integer(kind=jpim) :: jbegin_vder_EW = 0
integer(kind=jpim) :: jend_vder_EW = 0

logical :: ldump_values = .false.

integer, external :: ec_mpirank
logical :: luse_mpi = .true.

real(kind=jprb) :: zexwn, zeywn

!===================================================================================================

#include "setup_trans0.h"
#include "esetup_trans.h"
#include "einv_trans.h"
#include "edir_trans.h"
#include "etrans_inq.h"
#include "especnorm.h"
#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"
#include "ec_meminfo.intfb.h"

!===================================================================================================

luse_mpi = detect_mpirun()

! Setup
call get_command_line_arguments(nlon, nlat, nsmax, nmsmax, iters, nfld, nlev, lvordiv, lscders, luvders, &
  & nproma, verbosity, ldump_values, lprint_norms, lmeminfo, nprgpns, nprgpew, nprtrv, nprtrw, ncheck)
! derived defaults
if ( nsmax == 0 ) nsmax = nlat/2-1
if ( nmsmax == 0 ) nmsmax = nlon/2-1
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
endif
nthread = oml_max_threads()

call dr_hook_init()

!===================================================================================================

if( lstats ) call gstats(0,0)
ztinit = omp_get_wtime()

! only output to stdout on pe 1
!if (nproc > 1) then
  !if (myproc /= 1) then
    !open(unit=nout, file='output_'//char(myproc/10+48)//char(myproc+48)//'.dat')
  !endif
!endif

if (ldetailed_stats) then
  lstats_omp    = .true.
  lstats_comms  = .true.
  lstats_mpl    = .true.
  lstatscpu     = .true.
  nprnt_stats   = nproc
  lstats_mem   = .true.
  lstats_alloc = .true.
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
if (nproc == 0) nproc = 1
if ( nprgpew == 0 .and. nprgpns == 0 ) then
  isqr = int(sqrt(real(nproc,jprb)))
  do ja = isqr, nproc
    ib = nproc/ja
    if (ja*ib == nproc) then
    nprgpns = max(ja,ib)
    nprgpew = min(ja,ib)
    exit
    endif
  enddo
elseif (nprgpns == 0 ) then
  nprgpns=nproc/nprgpew
elseif (nprgpew == 0 ) then
  nprgpew=nproc/nprgpns
endif
if (nprgpns*nprgpew /= nproc) call abor1('transform_test:nprgpns*nprgpew /= nproc')

! From sumpini, although this should be specified in namelist
if (nspecresmin == 0) nspecresmin = nproc

! Compute nprtrv and nprtrw if not provided on the command line
if (nprtrv ==0 .and. nprtrw == 0 ) then
  nprtrv=nprgpew
  nprtrw=nprgpns
elseif (nprtrv == 0 ) then
  nprtrv=nproc/nprtrw
elseif (nprtrw == 0 ) then
  nprtrw=nproc/nprtrv
endif
if (nprtrv*nprtrw /= nproc) call abor1('transform_test:nprtrv*nprtrw /= nproc')

mysetv=mod(myproc-1,nprtrv)+1

! Determine number of local levels for zonal and meridional fourier calculations
! based on the values of nflevg and nprtrv
allocate(numll(nprtrv))
numll=nflevg/nprtrv
numll(1:modulo(nflevg,nprtrv))=numll(1:modulo(nflevg,nprtrv))+1
ivsetsc(1)=min(nflevg+1, nprtrv)
nflevl = numll(mysetv)

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

if( lstats ) call gstats(1, 0)
call setup_trans0(kout=nout, kerr=nerr, kprintlev=merge(2, 0, verbosity == 1),                &
  &               kmax_resol=nmax_resol, kpromatr=0, kprgpns=nprgpns, kprgpew=nprgpew, &
  &               kprtrw=nprtrw, ldsync_trans=lsync_trans,               &
  &               ldalloperm=.true., ldmpoff=.not.luse_mpi)
  if( lstats ) call gstats(1, 1)

  if( lstats ) call gstats(2, 0)
zexwn=1._jprb  ! 2*pi/(nx*dx): spectral resolution
zeywn=1._jprb  ! 2*pi/(ny*dy)
nloen=nlon
call esetup_trans(ksmax=nsmax, kmsmax=nmsmax, kdgl=nlat, kdgux=nlat, kloen=nloen, ldsplit=.true.,          &
  &                 ldusefftw=lfftw,pexwn=zexwn,peywn=zeywn)

  if( lstats ) call gstats(2, 1)

call etrans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg)

if (nproma == 0) then ! no blocking (default when not specified)
  nproma = ngptot
endif

! Calculate number of NPROMA blocks
ngpblks = (ngptot - 1)/nproma+1

!===================================================================================================
! Print information before starting
!===================================================================================================

! Print configuration details
if (verbosity >= 0) then
  write(nout,'(" ")')
  write(nout,'(a)')'======= Start of runtime parameters ======='
  write(nout,'(" ")')
  write(nout,'("nlon      ",i0)') nlon
  write(nout,'("nlat      ",i0)') nlat
  write(nout,'("nsmax     ",i0)') nsmax
  write(nout,'("nmsmax    ",i0)') nmsmax
  write(nout,'("nproc     ",i0)') nproc
  write(nout,'("nthread   ",i0)') nthread
  write(nout,'("nprgpns   ",i0)') nprgpns
  write(nout,'("nprgpew   ",i0)') nprgpew
  write(nout,'("nprtrw    ",i0)') nprtrw
  write(nout,'("nprtrv    ",i0)') nprtrv
  write(nout,'("ngptot    ",i0)') ngptot
  write(nout,'("ngptotg   ",i0)') ngptotg
  write(nout,'("nfld      ",i0)') nfld
  write(nout,'("nlev      ",i0)') nlev
  write(nout,'("nflevl    ",i0)') nflevl
  write(nout,'("nproma    ",i0)') nproma
  write(nout,'("ngpblks   ",i0)') ngpblks
  write(nout,'("nspec2    ",i0)') nspec2
  write(nout,'("nspec2g   ",i0)') nspec2g
  write(nout,'("lvordiv   ",l)') lvordiv
  write(nout,'("lscders   ",l)') lscders
  write(nout,'("luvders   ",l)') luvders
  write(nout,'(" ")')
  write(nout,'(a)') '======= End of runtime parameters ======='
  write(nout,'(" ")')
end if

!===================================================================================================
! Allocate and Initialize spectral arrays
!===================================================================================================

! Allocate spectral arrays
! Try to mimick IFS layout as much as possible
nullify(zspvor)
nullify(zspdiv)
nullify(zspsc3a)
allocate(sp3d(nflevl,nspec2,2+nfld))
allocate(zspsc2(1,nspec2))
allocate(zmeanu(nflevl),zmeanv(nflevl))
zmeanu(:)=0._jprb
zmeanv(:)=0._jprb

call initialize_spectral_arrays(nsmax, nmsmax, zspsc2, sp3d)

! Point convenience variables to storage variable sp3d
zspvor  => sp3d(:,:,1)
zspdiv  => sp3d(:,:,2)
zspsc3a => sp3d(:,:,3:3+(nfld-1))

!===================================================================================================
! Allocate gridpoint arrays
!===================================================================================================

allocate(ivset(nflevg))

! Compute spectral distribution
ilev = 0
do jb = 1, nprtrv
  do jlev=1, numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo

! Allocate grid-point arrays
if (lvordiv) then
  jbegin_uv = 1
  jend_uv = 2
endif
if (luvders) then
  jbegin_uder_EW  = jend_uv + 1
  jend_uder_EW    = jbegin_uder_EW + 1
  jbegin_vder_EW  = jend_uder_EW + 1
  jend_vder_EW    = jbegin_vder_EW + 1
else
  jbegin_uder_EW = jend_uv
  jend_uder_EW   = jend_uv
  jbegin_vder_EW = jend_uv
  jend_vder_EW   = jend_uv
endif

jbegin_sc = jend_vder_EW + 1
jend_sc   = jend_vder_EW + nfld

if (lscders) then
  ndimgmvs = 3
  jbegin_scder_NS = jend_sc + 1
  jend_scder_NS   = jend_sc + nfld
  jbegin_scder_EW = jend_scder_NS + 1
  jend_scder_EW   = jend_scder_NS + nfld
else
  ndimgmvs = 1
  jbegin_scder_NS = jend_sc
  jend_scder_NS   = jend_sc
  jbegin_scder_EW = jend_sc
  jend_scder_EW   = jend_sc
endif

ndimgmv = jend_scder_EW

!allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
!allocate(zgmvs(nproma,ndimgmvs,ngpblks))
!zgpuv => zgmv(:,:,1:jend_vder_EW,:)
!zgp3a => zgmv(:,:,jbegin_sc:jend_scder_EW,:)
!zgp2  => zgmvs(:,:,:)

! allocate separately since non-contiguous host-device transfers are not supported.
allocate(zgpuv(nproma,nflevg,jend_vder_EW,ngpblks))
allocate(zgp3a(nproma,nflevg,jend_scder_EW-jbegin_sc+1,ngpblks))
allocate(zgp2(nproma,ndimgmvs,ngpblks))

zgp2=0.
zgp3a=0.
zgpuv=0.

!===================================================================================================
! Allocate norm arrays
!===================================================================================================

if (lprint_norms .or. ncheck > 0) then
  allocate(znormsp(1))
  allocate(znormsp0(1))
  allocate(znormvor(nflevg))
  allocate(znormvor0(nflevg))
  allocate(znormdiv(nflevg))
  allocate(znormdiv0(nflevg))
  allocate(znormt(nflevg))
  allocate(znormt0(nflevg))

  call especnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor0, kvset=ivset(1:nflevg))
  call especnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv0, kvset=ivset(1:nflevg))
  call especnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt0,   kvset=ivset(1:nflevg))
  call especnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp0,  kvset=ivsetsc)
  
  if (verbosity >= 1 .and. myproc == 1) then
    do ifld = 1, nflevg
      write(nout,'("norm zspvor( ",i4,",:)   = ",f20.15)') ifld, znormvor0(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv0(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormt0(ifld)
    enddo
    do ifld = 1, 1
      write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15)') ifld, znormsp0(ifld)
    enddo
  endif
endif

!===================================================================================================
! Setup timers
!===================================================================================================

ztinit = (omp_get_wtime() - ztinit)

if (verbosity >= 0) then
  write(nout,'(" ")')
  write(nout,'(a,i6,a,f9.2,a)') "transform_test initialisation, on",nproc,&
                                & " tasks, took",ztinit," sec"
  write(nout,'(" ")')
endif

if (iters <= 0) call abor1('transform_test:iters <= 0')

allocate(ztstep(iters))
allocate(ztstep1(iters))
allocate(ztstep2(iters))

ztstepavg  = 0._jprd
ztstepmax  = 0._jprd
ztstepmin  = 9999999999999999._jprd
ztstepavg1 = 0._jprd
ztstepmax1 = 0._jprd
ztstepmin1 = 9999999999999999._jprd
ztstepavg2 = 0._jprd
ztstepmax2 = 0._jprd
ztstepmin2 = 9999999999999999._jprd

!=================================================================================================
! Dump the values to disk, for debugging only
!=================================================================================================

if (ldump_values) then
	! dump a field to a binary file
	call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspsc2(1,:),ivsetsc(1:1), 'S', noutdump)
	if (lvordiv) then
	  call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspdiv(1,:),ivset(1:1), 'D', noutdump)
	  call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspvor(1,:),ivset(1:1), 'V', noutdump)
	endif
    call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspsc3a(1,:,1),ivset(1:1), 'T', noutdump)
endif

write(nout,'(a)') '======= Start of spectral transforms  ======='
write(nout,'(" ")')

ztloop = omp_get_wtime()

!===================================================================================================
! Do spectral transform loop
!===================================================================================================

do jstep = 1, iters
  if( lstats ) call gstats(3,0)
  ztstep(jstep) = omp_get_wtime()

  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

  ztstep1(jstep) = omp_get_wtime()
  if( lstats ) call gstats(4,0)
  if (lvordiv) then

    call einv_trans(kresol=1, kproma=nproma, &
       & pspsc2=zspsc2,                     & ! spectral surface pressure
       & pspvor=zspvor,                     & ! spectral vorticity
       & pspdiv=zspdiv,                     & ! spectral divergence
       & pspsc3a=zspsc3a,                   & ! spectral scalars
       & ldscders=lscders,                  &
       & ldvorgp=.false.,                   & ! no gridpoint vorticity
       & lddivgp=.false.,                   & ! no gridpoint divergence
       & lduvder=luvders,                   &
       & kvsetuv=ivset,                     &
       & kvsetsc2=ivsetsc,                  &
       & kvsetsc3a=ivset,                   &
       & pgp2=zgp2,                         &
       & pgpuv=zgpuv,                       &
       & pgp3a=zgp3a,                       &
     & pmeanu=zmeanu,                     &
     & pmeanv=zmeanv)

  else

    call einv_trans(kresol=1, kproma=nproma, &
       & pspsc2=zspsc2,                     & ! spectral surface pressure
       & pspsc3a=zspsc3a,                   & ! spectral scalars
       & ldscders=lscders,                  & ! scalar derivatives
       & kvsetsc2=ivsetsc,                  &
       & kvsetsc3a=ivset,                   &
       & pgp2=zgp2,                         &
       & pgp3a=zgp3a)

  endif
  
  if( lstats ) call gstats(4,1)

  ztstep1(jstep) = (omp_get_wtime() - ztstep1(jstep))

  !=================================================================================================
  ! While in grid point space, dump the values to disk, for debugging only
  !=================================================================================================

  if (ldump_values) then
    ! dump a field to a binary file
    call dump_gridpoint_field(jstep, myproc, nlat, nproma, ngpblks, zgp2(:,1,:),         'S', noutdump)
    if (lvordiv) then
        call dump_gridpoint_field(jstep, myproc, nlat, nproma, ngpblks, zgpuv(:,nflevg,1,:), 'U', noutdump)
        call dump_gridpoint_field(jstep, myproc, nlat, nproma, ngpblks, zgpuv(:,nflevg,2,:), 'V', noutdump)
    endif
    call dump_gridpoint_field(jstep, myproc, nlat, nproma, ngpblks, zgp3a(:,nflevg,1,:), 'T', noutdump)
  endif
  
  !=================================================================================================
  ! Do direct transform
  !=================================================================================================

  ztstep2(jstep) = omp_get_wtime()

  if( lstats ) call gstats(5,0)
  

  if (lvordiv) then
    call edir_trans(kresol=1, kproma=nproma, &
      & pgp2=zgp2(:,1:1,:),                &
      & pgpuv=zgpuv(:,:,1:2,:),             &
      & pgp3a=zgp3a(:,:,1:nfld,:),          &
      & pspvor=zspvor,                      &
      & pspdiv=zspdiv,                      &
      & pspsc2=zspsc2,                      &
      & pspsc3a=zspsc3a,                    &
      & kvsetuv=ivset,                      &
      & kvsetsc2=ivsetsc,                   &
      & kvsetsc3a=ivset,                    &
    & pmeanu=zmeanu,                      &
    & pmeanv=zmeanv)
  else
  
    call edir_trans(kresol=1, kproma=nproma, &
      & pgp2=zgp2(:,1:1,:),                &
      & pgp3a=zgp3a(:,:,1:nfld,:),          &
      & pspsc2=zspsc2,                      &
      & pspsc3a=zspsc3a,                    &
      & kvsetsc2=ivsetsc,                   &
      & kvsetsc3a=ivset)
  endif
  if( lstats ) call gstats(5,1)
  ztstep2(jstep) = (omp_get_wtime() - ztstep2(jstep))

  !=================================================================================================
	! Dump the values to disk, for debugging only
	!=================================================================================================

	if (ldump_values) then
		! dump a field to a binary file
		call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspsc2(1,:),ivsetsc(1:1), 'S', noutdump)
		if (lvordiv) then
		  call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspdiv(1,:),ivset(1), 'D', noutdump)
		  call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspvor(1,:),ivset(1:1), 'V', noutdump)
		endif
		call dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, zspsc3a(1,:,1),ivset(1:1), 'T', noutdump)
	endif

  !=================================================================================================
  ! Calculate timings
  !=================================================================================================

  ztstep(jstep) = (omp_get_wtime() - ztstep(jstep))

  ztstepavg = ztstepavg + ztstep(jstep)
  ztstepmin = min(ztstep(jstep), ztstepmin)
  ztstepmax = max(ztstep(jstep), ztstepmax)

  ztstepavg1 = ztstepavg1 + ztstep1(jstep)
  ztstepmin1 = min(ztstep1(jstep), ztstepmin1)
  ztstepmax1 = max(ztstep1(jstep), ztstepmax1)

  ztstepavg2 = ztstepavg2 + ztstep2(jstep)
  ztstepmin2 = min(ztstep2(jstep), ztstepmin2)
  ztstepmax2 = max(ztstep2(jstep), ztstepmax2)

  !=================================================================================================
  ! Print norms
  !=================================================================================================

  if (lprint_norms) then
    if( lstats ) call gstats(6,0)
    call especnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp,  kvset=ivsetsc(1:1))
    call especnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor, kvset=ivset(1:nflevg))
    call especnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv, kvset=ivset(1:nflevg))
    call especnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt,   kvset=ivset(1:nflevg))
  
    if ( myproc == 1 ) then

      ! Surface pressure
      zmaxerr(:) = -999.0
      do ifld = 1, 1
        zerr(1) = abs(znormsp(ifld)/znormsp0(ifld) - 1.0_jprb)
        zmaxerr(1) = max(zmaxerr(1), zerr(1))
      enddo
      ! Divergence
      do ifld = 1, nflevg
        zerr(2) = abs(znormdiv(ifld)/znormdiv0(ifld) - 1.0_jprb)
        zmaxerr(2) = max(zmaxerr(2), zerr(2))
      enddo
      ! Vorticity
      do ifld = 1, nflevg
        zerr(3) = abs(znormvor(ifld)/znormvor0(ifld) - 1.0_jprb)
        zmaxerr(3) = max(zmaxerr(3),zerr(3))
      enddo
      ! Temperature
      do ifld = 1, nflevg
        zerr(4) = abs(znormt(ifld)/znormt0(ifld) - 1.0_jprb)
        zmaxerr(4) = max(zmaxerr(4), zerr(4))
      enddo
      write(nout,'("time step ",i6," took", f8.4," | zspvor max err="e10.3,&
                  & " | zspdiv max err="e10.3," | zspsc3a max err="e10.3," | zspsc2 max err="e10.3)') &
                  &  jstep, ztstep(jstep), zmaxerr(3), zmaxerr(2), zmaxerr(4), zmaxerr(1)
      if( lstats )call gstats(6,1)
    else
      write(nout,'("Time step ",i6," took", f8.4)') jstep, ztstep(jstep)
    endif

	endif

  if( lstats ) call gstats(3,1)

enddo

!===================================================================================================

ztloop = (omp_get_wtime() - ztloop)

write(nout,'(" ")')
write(nout,'(a)') '======= End of spectral transforms  ======='
write(nout,'(" ")')

if (lprint_norms .or. ncheck > 0) then
  call especnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor, kvset=ivset)
  call especnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv, kvset=ivset)
  call especnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt,   kvset=ivset)
  call especnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp,  kvset=ivsetsc)
  
  if ( myproc == 1 ) then

	  zmaxerr(:) = -999.0
	  do ifld = 1, nflevg
		zerr(3) = abs(real(znormvor(ifld),kind=jprd)/real(znormvor0(ifld),kind=jprd) - 1.0_jprd)
		zmaxerr(3) = max(zmaxerr(3), zerr(3))
		if (verbosity >= 1) then
		  write(nout,'("norm zspvor( ",i4,")     = ",f20.15,"        error = ",e10.3)') ifld, znormvor0(ifld), zerr(3)
		endif
	  enddo
	  do ifld = 1, nflevg
		zerr(2) = abs(real(znormdiv(ifld),kind=jprd)/real(znormdiv0(ifld),kind=jprd) - 1.0d0)
		zmaxerr(2) = max(zmaxerr(2),zerr(2))
		if (verbosity >= 1) then
		  write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15,"        error = ",e10.3)') ifld, znormdiv0(ifld), zerr(2)
		endif
	  enddo
	  do ifld = 1, nflevg
		zerr(4) = abs(real(znormt(ifld),kind=jprd)/real(znormt0(ifld),kind=jprd) - 1.0d0)
		zmaxerr(4) = max(zmaxerr(4), zerr(4))
		if (verbosity >= 1) then
		  write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15,"        error = ",e10.3)') ifld, znormt0(ifld), zerr(4)
		endif
	  enddo
	  do ifld = 1, 1
		zerr(1) = abs(real(znormsp(ifld),kind=jprd)/real(znormsp0(ifld),kind=jprd) - 1.0d0)
		zmaxerr(1) = max(zmaxerr(1), zerr(1))
		if (verbosity >= 1) then
		  write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15,"        error = ",e10.3)') ifld, znormsp0(ifld), zerr(1)
		endif
	  enddo

	  ! maximum error across all fields
	  zmaxerrg = max(max(zmaxerr(1),zmaxerr(2)), max(zmaxerr(2), zmaxerr(3)))

	  if (verbosity >= 1) write(nout,*)
	  write(nout,'("max error zspvor(1:nlev,:)    = ",e10.3)') zmaxerr(3)
	  write(nout,'("max error zspdiv(1:nlev,:)    = ",e10.3)') zmaxerr(2)
	  write(nout,'("max error zspsc3a(1:nlev,:,1) = ",e10.3)') zmaxerr(4)
	  write(nout,'("max error zspsc2(1:1,:)       = ",e10.3)') zmaxerr(1)
	  write(nout,*)
	  write(nout,'("max error combined =          = ",e10.3)') zmaxerrg
	  write(nout,*)

	  if (ncheck > 0) then
		! If the maximum spectral norm error across all fields is greater than 100 times the machine
		! epsilon, fail the test
		if (zmaxerrg > real(ncheck, jprb) * epsilon(1.0_jprb)) then
		  write(nout, '(a)') '*******************************'
		  write(nout, '(a)') 'Correctness test failed'
		  write(nout, '(a,1e7.2)') 'Maximum spectral norm error = ', zmaxerrg
		  write(nout, '(a,1e7.2)') 'Error tolerance = ', real(ncheck, jprb) * epsilon(1.0_jprb)
		  write(nout, '(a)') '*******************************'
		  error stop
		endif
	  endif
  endif
endif

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

call sort(ztstep,iters)
ztstepmed = ztstep(iters/2)

ztstepavg1 = (ztstepavg1/real(nproc,jprb))/real(iters,jprd)
ztstep1(:) = ztstep1(:)/real(nproc,jprd)

call sort(ztstep1, iters)
ztstepmed1 = ztstep1(iters/2)

ztstepavg2 = (ztstepavg2/real(nproc,jprb))/real(iters,jprd)
ztstep2(:) = ztstep2(:)/real(nproc,jprd)

call sort(ztstep2,iters)
ztstepmed2 = ztstep2(iters/2)


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
      call mpl_recv(istack, ksource=nprcids(i), ktag=i, cdstring='transform_test:')
      print '(i4,11x,i10)', i, istack
    enddo
  else
    call mpl_send(istack, kdest=nprcids(1), ktag=myproc, cdstring='transform_test:')
  endif
endif

!===================================================================================================
! Cleanup
!===================================================================================================

! TODO: many more arrays to deallocate

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

subroutine get_command_line_arguments(nlon, nlat, nsmax, nmsmax, &
 &                                    iters, nfld, nlev, lvordiv, lscders, luvders, &
  &                                   nproma, verbosity, ldump_values, lprint_norms, &
  &                                   lmeminfo, nprgpns, nprgpew, nprtrv, nprtrw, ncheck)

  integer, intent(inout) :: nlon            ! Zonal dimension
  integer, intent(inout) :: nlat            ! Meridional dimension
  integer, intent(inout) :: nsmax           ! Meridional truncation
  integer, intent(inout) :: nmsmax          ! Zonal trunciation
  integer, intent(inout) :: iters           ! Number of iterations for transform test
  integer, intent(inout) :: nfld            ! Number of scalar fields
  integer, intent(inout) :: nlev            ! Number of vertical levels
  logical, intent(inout) :: lvordiv         ! Also transform vorticity/divergence
  logical, intent(inout) :: lscders         ! Compute scalar derivatives
  logical, intent(inout) :: luvders         ! Compute uv East-West derivatives
  integer, intent(inout) :: nproma          ! NPROMA
  integer, intent(inout) :: verbosity       ! Level of verbosity
  logical, intent(inout) :: ldump_values    ! Dump values of grid point fields for debugging
  logical, intent(inout) :: lprint_norms    ! Calculate and print spectral norms of fields
  logical, intent(inout) :: lmeminfo        ! Show information from FIAT ec_meminfo routine at the
                                            ! end
  integer, intent(inout) :: nprgpns         ! Size of NS set (gridpoint decomposition)
  integer, intent(inout) :: nprgpew         ! Size of EW set (gridpoint decomposition)
  integer, intent(inout) :: nprtrv          ! Size of V set (spectral decomposition)
  integer, intent(inout) :: nprtrw          ! Size of W set (spectral decomposition)
  integer, intent(inout) :: ncheck          ! The multiplier of the machine epsilon used as a
                                            ! tolerance for correctness checking

  character(len=128) :: carg          ! Storage variable for command line arguments
  integer            :: iarg = 1      ! Argument index
  integer            :: stat          ! For storing success status of string->integer conversion
  integer            :: myproc

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
      ! Parse spectral truncation argument
      case('--nlon'); nlon = get_int_value('--nlon', iarg)
      case('--nlat'); nlat = get_int_value('--nlat', iarg)
      case('--nsmax'); nsmax = get_int_value('--nsmax', iarg)
      case('--nmsmax'); nmsmax = get_int_value('--nmsmax', iarg)
      case('-f', '--nfld'); nfld = get_int_value('-f', iarg)
      case('-l', '--nlev'); nlev = get_int_value('-l', iarg)
      case('--vordiv'); lvordiv = .True.
      case('--scders'); lscders = .True.
      case('--uvders'); luvders = .True.
      case('--nproma'); nproma = get_int_value('--nproma', iarg)
      case('--dump-values'); ldump_values = .true.
      case('--norms'); lprint_norms = .true.
      case('--meminfo'); lmeminfo = .true.
      case('--nprgpns'); nprgpns = get_int_value('--nprgpns', iarg)
      case('--nprgpew'); nprgpew = get_int_value('--nprgpew', iarg)
      case('--nprtrv'); nprtrv = get_int_value('--nprtrv', iarg)
      case('--nprtrw'); nprtrw = get_int_value('--nprtrw', iarg)
      case('-c', '--check'); ncheck = get_int_value('-c', iarg)
      case default
        call parsing_failed("Unrecognised argument: " // trim(carg))

    end select
    iarg = iarg + 1
  end do

  if (.not. lvordiv) then
    luvders = .false.
  endif

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

subroutine print_help(unit)

  integer, optional :: unit
  integer :: nout = 6
  if (present(unit)) then
    nout = unit
  endif

  write(nout, "(a)") ""

  if (jprb == jprd) then
    write(nout, "(a)") "NAME    ectrans-lam-benchmark-dp"
  else
    write(nout, "(a)") "NAME    ectrans-lam-benchmark-sp"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "DESCRIPTION"
  write(nout, "(a)") "        This program tests ecTrans-lam by transforming fields back and forth&
    & between spectral "
  if (jprb == jprd) then
    write(nout, "(a)") "        space and grid-point space (double-precision version)"
  else
    write(nout, "(a)") "        space and grid-point space (single-precision version)"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "USAGE"
  if (jprb == jprd) then
    write(nout, "(a)") "        ectrans-lam-benchmark-dp [options]"
  else
    write(nout, "(a)") "        ectrans-lam-benchmark-sp [options]"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "OPTIONS"
  write(nout, "(a)") "    -h, --help          Print this message"
  write(nout, "(a)") "    -v                  Run with verbose output"
  write(nout, "(a)") "    --nlon NLON         Number of gridpoints in zonal direction (default = 128)"
  write(nout, "(a)") "    --nlat NLAT         Number of gridpoints in meridional direction (default = 128)"
  write(nout, "(a)") "    --nsmax NSMAX       Spectral truncation in meridional direction (default = NLAT/2-1)"
  write(nout, "(a)") "    --nmsmax NMSMAX     Spectral truncation in zonal direction (default = NLON/2-1)"
  write(nout, "(a)") "    -n, --niter NITER   Run for this many inverse/direct transform&
    & iterations (default = 10)"
  write(nout, "(a)") "    -f, --nfld NFLD     Number of scalar fields (default = 1)"
  write(nout, "(a)") "    -l, --nlev NLEV     Number of vertical levels (default = 1)"
  write(nout, "(a)") "    --vordiv            Also transform vorticity-divergence to wind"
  write(nout, "(a)") "    --scders            Compute scalar derivatives (default off)"
  write(nout, "(a)") "    --uvders            Compute uv East-West derivatives (default off). Only&
    & when also --vordiv is given"
  write(nout, "(a)") "    --nproma NPROMA     Run with NPROMA (default no blocking: NPROMA=ngptot)"
  write(nout, "(a)") "    --norms             Calculate and print spectral norms of transformed&
    & fields"
  write(nout, "(a)") "                        The computation of spectral norms will skew overall&
    & timings"
  write(nout, "(a)") "    --meminfo           Show diagnostic information from FIAT's ec_meminfo&
    & subroutine on memory usage, thread-binding etc."
  write(nout, "(a)") "    --nprgpew           Size of East-West set in gridpoint decomposition"
  write(nout, "(a)") "    --nprgpns           Size of North-South set in gridpoint decomposition"
  write(nout, "(a)") "    --nprtrv            Size of Vertical set in spectral decomposition"
  write(nout, "(a)") "    --nprtrw            Size of Wave set in spectral decomposition"
  write(nout, "(a)") "    -c, --check VALUE   The multiplier of the machine epsilon used as a&
   & tolerance for correctness checking"
  write(nout, "(a)") ""
  write(nout, "(a)") "DEBUGGING"
  write(nout, "(a)") "    --dump-values       Output gridpoint fields in unformatted binary file"
  write(nout, "(a)") ""

end subroutine print_help

!===================================================================================================

subroutine initialize_spectral_arrays(nsmax, nmsmax, zsp, sp3d)

  integer,         intent(in)    :: nsmax     ! Spectral truncation in meridional direction
  integer,         intent(in)    :: nmsmax    ! Spectral truncation in zonal direction
  real(kind=jprb), intent(inout) :: zsp(:,:)    ! Surface pressure
  real(kind=jprb), intent(inout) :: sp3d(:,:,:) ! 3D fields

  integer(kind=jpim) :: nflevl
  integer(kind=jpim) :: nfield

  integer :: i, j

  nflevl = size(sp3d, 1)
  nfield = size(sp3d, 3)

  ! First initialize surface pressure
  call initialize_2d_spectral_field(nsmax, nmsmax, zsp(1,:))

  ! Then initialize all of the 3D fields
  do i = 1, nflevl
    do j = 1, nfield
      call initialize_2d_spectral_field(nsmax, nmsmax, sp3d(i,:,j))
    end do
  end do

end subroutine initialize_spectral_arrays

!===================================================================================================

subroutine initialize_2d_spectral_field(nsmax, nmsmax, field)

  integer,         intent(in)    :: nsmax     ! Spectral truncation in meridional direction
  integer,         intent(in)    :: nmsmax    ! Spectral truncation in zonal direction
  real(kind=jprb), intent(inout) :: field(:)  ! Field to initialize

  integer :: ispec, kspec2
  integer, allocatable :: my_km(:), my_kn(:)

  ! Choose a harmonic to initialize arrays
  integer :: m_num = 1 ! Zonal wavenumber
  integer :: n_num = 0 ! Meridional wavenumber
  
  ! Type of initialization: (single) 'harmonic' or (random) 'spectrum'
  character(len=32) :: init_type='harmonic'    

  ! First initialise all spectral coefficients to zero
  field(:) = 0.0
  
  ! make sure wavenumbers are within truncation
  if ( m_num>nmsmax .or. n_num > nsmax .or. &
     & ( nsmax>0 .and. nmsmax>0 .and. ( (m_num/real(nmsmax))**2+(n_num/real(nsmax))**2 ) > 1.) ) then
    write (nerr,*)
    write (nerr,*) 'WARNING: INITIAL WAVENUMBERS OUTSIDE OF TRUNCATION! '
    write (nerr,*) '  m_num = ',m_num,'; nmsmax = ',nmsmax,'; n_num = ',n_num,'; nsmax = ',nsmax,&
     & '; ellips check: ',(m_num/real(nmsmax))**2+(n_num/real(nsmax))**2
    write (nerr,*) '  using (kx=',NMSMAX/2,', ky=', NSMAX/2,') instead'
    write (nerr,*)
    m_num=nmsmax/2
    n_num=nsmax/2
  endif
  
  ! Get wavenumbers this rank is responsible for
  call etrans_inq(kspec2=kspec2)
  allocate(my_kn(kspec2),my_km(kspec2))
  call etrans_inq(knvalue=my_kn,kmvalue=my_km)

  ! If rank is responsible for the chosen zonal wavenumber...
  if ( init_type == 'harmonic' ) then
    do ispec=1,nspec2,4
      if ( my_kn(ispec)== n_num .and. my_km(ispec) == m_num ) then
        field(ispec)=1.0 ! cos*cos
        !field(ispec+1)=1.0 ! cos*sin
        !field(ispec+2)=1.0 ! sin*cos
        !field(ispec+3)=1.0 ! sin*sin
      end if
    enddo
  endif

  ! random power spectrum
  if ( init_type == 'spectrum' ) then
    call random_number(field)
    field=2*field-1.  ! center around zero
    ! set some components to zero because they are unphysical
    do ispec=1,nspec2,4
      if ( my_kn(ispec)== 0 .and. my_km(ispec) == 0 ) field(ispec:ispec+3)=0. ! remove mean value for vorticity and divergence
      if ( my_kn(ispec)== 0 ) field(ispec+1:ispec+3:2)=0. ! remove sine component on zero-wavenumber
      if ( my_kn(ispec)== nmsmax ) field(ispec+1:ispec+3:2)=0. ! remove sine component on last-wavenumber
      if ( my_km(ispec)== 0 ) field(ispec+2:ispec+3)=0. ! remove sine component on zero-wavenumber
      if ( my_km(ispec)== nsmax ) field(ispec+2:ispec+3)=0. ! remove sine component on last-wavenumber
    enddo
    
    ! scale according to wavenumber**2
    do ispec=1,nspec2
      field(ispec)=field(ispec)/(0.01+(my_kn(ispec)/real(nsmax))**2+(my_km(ispec)/real(nmsmax))**2)
    enddo
  endif
  
end subroutine initialize_2d_spectral_field

!===================================================================================================

subroutine dump_gridpoint_field(jstep, myproc, nlat, nproma, ngpblks, fld, fldchar, noutdump)

  ! Dump a 2d gridpoint field to screen or a binary file.

  integer(kind=jpim), intent(in) :: jstep ! Time step, used for naming file
  integer(kind=jpim), intent(in) :: myproc ! MPI rank, used for naming file
  integer(kind=jpim), intent(in) :: nlat ! Number of latitudes
  integer(kind=jpim), intent(in) :: nproma ! Size of nproma
  integer(kind=jpim), intent(in) :: ngpblks ! Number of nproma blocks
  real(kind=jprb)   , intent(in) :: fld(nproma,1,ngpblks) ! 2D field
  character         , intent(in) :: fldchar ! Single character field identifier
  integer(kind=jpim), intent(in) :: noutdump ! Unit number for output file
  
  integer(kind=jpim) :: kgptotg      ! global number of gridpoints
  real(kind=jprb), allocatable :: fldg(:,:)  ! global field
  integer(kind=jpim) :: kfgathg=1    ! number of fields to gather
  integer(kind=jpim) :: kto(1)=(/1/) ! processor where to gather
  character(len=14)  :: filename = "x.xxx.xxx.grid"
  character(len=13) :: frmt='(4X,xxxxF8.2)'

#include "etrans_inq.h"
#include "egath_grid.h"
  
  call etrans_inq(kgptotg=kgptotg)

  if ( myproc == 1 ) allocate(fldg(kgptotg,1))

  call egath_grid(pgpg=fldg,kproma=nproma,kfgathg=kfgathg,kto=kto,pgp=fld)
  
  if ( myproc == 1 ) then
 
    ! write to file
    write(filename(1:1),'(a1)') fldchar
    write(filename(3:5),'(i3.3)') jstep
#ifdef ACCGPU
    write(filename(7:9),'(a3)') 'gpu'
#else
    write(filename(7:9),'(a3)') 'cpu'
#endif
    open(noutdump, file=filename, form="unformatted", access="stream")
    write(noutdump) kgptotg/nlat,nlat ! dimensions
	write(noutdump) fldg ! data
    close(noutdump)
    
    ! write to screen
    write(frmt(5:8),'(i4.4)') kgptotg/nlat
    write (*,*) fldchar,' at iteration ',jstep,':'
    write (*,frmt) fldg
    call flush(6)
	
    deallocate(fldg)
  
  endif
  

end subroutine dump_gridpoint_field

!===================================================================================================

subroutine dump_spectral_field(jstep, myproc, nspec2, nsmax, nmsmax, fld, kvset, fldchar, noutdump)

  ! Dump a 2d spectral field to screen or a binary file.

  integer(kind=jpim), intent(in) :: jstep ! Time step, used for naming file
  integer(kind=jpim), intent(in) :: myproc ! MPI rank, used for naming file
  integer(kind=jpim), intent(in) :: nspec2 ! Size of nspec2 (number of waves on this proc in M-space)
  integer(kind=jpim), intent(in) :: nsmax
  integer(kind=jpim), intent(in) :: nmsmax
  real(kind=jprb)   , intent(in) :: fld(1,nspec2) ! 2D field
  integer(kind=jpim), intent(in) :: kvset(1)   ! B-set on which the field resides
  character         , intent(in) :: fldchar ! Single character field identifier
  integer(kind=jpim), intent(in) :: noutdump ! Unit number for output file
  
  integer(kind=jpim) :: nspec2g              ! global number of gridpoints
  real(kind=jprb), allocatable :: fldg(:,:)  ! global field (nspec2g)
  integer(kind=jpim) :: kfgathg=1    ! number of fields to gather
  integer(kind=jpim) :: kto(1)=(/1/) ! processor where to gather
  character(len=14)   :: filename = "x.xxx.xxx.spec"
  character(len=13)  :: frmt='(4X,xxxxF8.2)' ! for printing to screen
  integer(kind=jpim) :: knse(0:nmsmax),kmse(0:nsmax) ! elliptic truncation
  real(kind=jprb)    :: fld2g(0:2*nmsmax+1,0:2*nsmax+1) ! 2D representation of spectral field
  integer(kind=jpim) :: jj, jms, jns
  
#include "etrans_inq.h"
#include "egath_spec.h"
  
  if ( myproc == 1 ) then
    call etrans_inq(kspec2g=nspec2g)
    allocate(fldg(1,nspec2g))
	call ellips(nsmax,nmsmax,knse,kmse)
  endif

  call egath_spec(PSPECG=fldg,kfgathg=kfgathg,kto=kto,kvset=kvset,PSPEC=fld)
  
  if ( myproc == 1 ) then

	fld2g=0.
	jj=1
	do jms=0,nmsmax
	  do jns=0,knse(jms)
		fld2g(2*jms+0,2*jns+0)=fldg(1,jj)
		fld2g(2*jms+0,2*jns+1)=fldg(1,jj+1)
		fld2g(2*jms+1,2*jns+0)=fldg(1,jj+2)
		fld2g(2*jms+1,2*jns+1)=fldg(1,jj+3)
		jj=jj+4
	  enddo
	enddo
	
    ! write to binary file
    write(filename(1:1),'(a1)') fldchar
    write(filename(3:5),'(i3.3)') jstep
#ifdef ACCGPU
    write(filename(7:9),'(a3)') 'gpu'
#else
    write(filename(7:9),'(a3)') 'cpu'
#endif
    open(noutdump, file=filename, form="unformatted", access="stream")
	write(noutdump) 2*nmsmax+2,2*nsmax+2  ! dimensions
    write(noutdump) fld2g             ! data
    close(noutdump)
    
    ! write to screen
    write(frmt(5:8),'(i4.4)') 2*(nmsmax+1)
    write (*,*) fldchar,' at iteration ',jstep,':'
    write (*,frmt) fld2g
    call flush(6)
	
    deallocate(fldg)
  
  endif
  

end subroutine dump_spectral_field

!===================================================================================================

function detect_mpirun() result(lmpi_required)
  logical :: lmpi_required
  integer :: ilen
  integer, parameter :: nvars = 5
  character(len=32), dimension(nvars) :: cmpirun_detect
  character(len=4) :: clenv_dr_hook_assert_mpi_initialized
  integer :: ivar

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  cmpirun_detect(1) = 'OMPI_COMM_WORLD_SIZE'  ! openmpi
  cmpirun_detect(2) = 'ALPS_APP_PE'           ! cray pe
  cmpirun_detect(3) = 'PMI_SIZE'              ! intel
  cmpirun_detect(4) = 'SLURM_NTASKS'          ! slurm
  cmpirun_detect(5) = 'ECTRANS_USE_MPI'       ! forced

  lmpi_required = .false.
  do ivar = 1, nvars
    call get_environment_variable(name=trim(cmpirun_detect(ivar)), length=ilen)
    if (ilen > 0) then
      lmpi_required = .true.
      exit ! break
    endif
  enddo
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

end program ectrans_lam_benchmark

!===================================================================================================