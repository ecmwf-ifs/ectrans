! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program ectrans_benchmark

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

use parkind1,              only : jpim, jprb, jprd
use mpl_module,            only : mpl_allreduce, mpl_broadcast

use ectrans_device_mod,    only : ectrans_device_is_host
use ectrans_program_mod,   only : ectrans_program_init, ectrans_program_end, ectrans_print_memory_usage, &
                                & nthread, nproc, myproc
use ectrans_log_mod,       only : nout, nerr
use ectrans_gstats_mod,    only : ectrans_gstats_enable
use ectrans_grids_mod,     only : parse_gaussian_grid, cubic_octahedral_gaussian_grid
use ectrans_mpi_mod,       only : ectrans_mpi_enabled
use ectrans_checksum_mod,  only : ectrans_checksum_file_writer
use ectrans_timer_mod,     only : ectrans_timer, ectrans_timings
use ectrans_decomposition_mod,    only : ectrans_spectral_decomposition, ectrans_gridpoint_decomposition, &
                                       & ectrans_make_spectral_distribution
use ectrans_memory_mod,    only : allocator
use ectrans_command_line_parser_mod, only : ectrans_command_line_parser

#if USE_FIELD_API
USE ectrans_field_api_helper, only : wrapped_fields, fields_lists, &
                                   & wrap_benchmark_fields, wrap_benchmark_fields_zgp, &
                                   & create_fields_lists, &
                                   & delete_wrapped_fields,delete_fields_lists, &
                                   & synchost_rdonly_wrapped_fields
#endif

implicit none

! Default parameters
integer(kind=jpim) :: iters   = 10  ! Number of iterations for transform test
integer(kind=jpim) :: nfld    = 1   ! Number of 3D scalar fields
integer(kind=jpim) :: nlev    = 1   ! Number of vertical levels
integer(kind=jpim) :: iters_warmup = 3 ! Number of warm up steps (for which timing statistics should be ignored)

integer(kind=jpim) :: nflevg  ! Total number of vertical levels
integer(kind=jpim) :: nflevl  ! Number of vertical levels in the V set

integer(kind=jpim) :: nspec2  ! Number of spectral coefficients (real and imaginary)
integer(kind=jpim) :: ngptot  ! Total number of grid points on this task
integer(kind=jpim) :: ngptotg ! Total number of grid points across all tasks
integer(kind=jpim) :: mysetv  ! V set index of this task (1-based)

integer(kind=jpim) :: ifld
integer(kind=jpim) :: jb
integer(kind=jpim) :: nspec2g

integer :: jstep

type(ectrans_timer) :: timer_loop
type(ectrans_timer) :: timer_tstep
type(ectrans_timer) :: timer_invtrans
type(ectrans_timer) :: timer_dirtrans
type(ectrans_timer) :: timer_init
type(ectrans_timer) :: timer_norms
type(ectrans_timer) :: timer_setup_trans0
type(ectrans_timer) :: timer_setup_trans
type(ectrans_timer) :: timer_dump_checksums
type(ectrans_timer) :: timer_dump_values
type(ectrans_timings) :: timings_loop
type(ectrans_timings) :: timings_tstep
type(ectrans_timings) :: timings_invtrans
type(ectrans_timings) :: timings_dirtrans

! Arrays to store norm computations
real(kind=jprb), allocatable :: znormvor(:), znormvor1(:), znormdiv(:), znormdiv1(:)
real(kind=jprb), allocatable :: znormscalar(:), znormscalar1(:)
real(kind=jprb), allocatable :: znormsc3a(:), znormsc3a1(:), znormsc2(:), znormsc21(:)

! Error checking arrays
real(kind=jprd) :: zmaxerr(4)
real(kind=jprd) :: zmaxerrg

! Spectral space data structures
real(kind=jprb), pointer :: zspvor(:,:)
real(kind=jprb), pointer :: zspdiv(:,:)
real(kind=jprb), pointer :: zspscalar(:,:)
real(kind=jprb), pointer :: zspsc3a(:,:,:)
real(kind=jprb), pointer :: zspsc2(:,:)

! Grid-point space data structures
real(kind=jprb), pointer :: zgp(:,:,:)
real(kind=jprb), pointer :: zgpuv(:,:,:,:)
real(kind=jprb), pointer :: zgp3a(:,:,:,:)
real(kind=jprb), pointer :: zgp2(:,:,:)

logical :: lfield_api = .false. ! Whether to use field API or not

! setup_trans options
integer(kind=jpim) :: nsmax   = 79  ! Spectral truncation
integer(kind=jpim) :: ndgl    ! Number of latitudes
integer(kind=jpim), allocatable :: nloen(:) ! Number of points on each latitude
logical :: luserpnm = .false. ! Use Belusov algorithm to compute RPNM array instead of per m
logical :: luseflt = .false. ! Use fast legendre transforms
integer(kind=jpim) :: nopt_mem_tr = 0
integer(kind=jpim) :: nprtrv = 0 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
integer(kind=jpim) :: nprgpns = 0 ! Gridpoint decomp
integer(kind=jpim) :: nprgpew = 0 ! Gridpoint decomp

! Extra inv_trans options
logical :: lvordiv = .false. ! Compute vorticity and divergence in grid point space
logical :: lscders = .false. ! Compute derivatives of scalar (North-South and East-West) in grid
                             ! point space
logical :: luvder = .false. ! Compute East-West derivatives of U and V wind in grid point space

logical :: lprint_norms = .false. ! Calculate and print spectral norms
logical :: lmeminfo = .false. ! Show information from FIAT routine ec_meminfo at the end

! The multiplier of the machine epsilon used as a tolerance for correctness checking
! ncheck = 0 (the default) means that correctness checking is disabled
integer(kind=jpim) :: ncheck = 0

! Verbosity level (0 or 1)
integer :: verbosity = 0

integer(kind=jpim), allocatable :: numll(:), ivset(:), ivsetsc(:)
integer(kind=jpim) :: ivsetsc2(1)


! sumpini
integer(kind=jpim) :: isqr
logical :: lsync_trans = .true. ! Activate barrier sync
logical :: leq_regions = .true. ! Eq regions flag

integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: npromatr = 0
integer(kind=jpim) :: ngpblks
integer(kind=jpim) :: iend

! locals
integer(kind=jpim) :: ilev, jlev, ierr

#if USE_FIELD_API
type(wrapped_fields) :: ywflds
type(fields_lists) :: ylf
#endif

logical :: ldump_values = .false.
logical :: lpinning = .false.
logical :: ldump_checksums = .false.
logical :: lalloperm = .true.

character(len=16)   :: cgrid = ''
character(len=128)  :: cchecksums_path = ''

integer :: icall_mode = 2
integer :: inum_wind_fields, inum_sc_3d_fields, inum_sc_2d_fields, itotal_fields
integer :: ipgp_start, ipgp_end, ipgpuv_start, ipgpuv_end, islice
logical :: lloop_timer_active = .false.

real(kind=jprb), allocatable :: global_field(:,:)

type(ectrans_checksum_file_writer) :: invtrans_checksums_file_writer
type(ectrans_checksum_file_writer) :: dirtrans_checksums_file_writer

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#if USE_FIELD_API
#include "inv_trans_field_api.h"
#include "dir_trans_field_api.h"
#endif
#include "trans_inq.h"
#include "gath_grid.h"
#include "specnorm.h"
#include "abor1.intfb.h"
#include "ec_meminfo.intfb.h"

!===================================================================================================

lpinning = .not. ectrans_device_is_host()

! Setup
call get_command_line_arguments(nsmax, cgrid, iters, iters_warmup, nfld, nlev, lvordiv, lscders, &
  &                             luvder, luseflt, nopt_mem_tr, nproma, npromatr, verbosity, &
  &                             ldump_values, lprint_norms, lmeminfo, nprtrv, nprtrw, ncheck, &
  &                             lpinning, lfield_api, icall_mode, ldump_checksums, cchecksums_path, lalloperm)
if (cgrid == '') cgrid = cubic_octahedral_gaussian_grid(nsmax)
call parse_gaussian_grid(cgrid, ndgl, nloen)
nflevg = nlev

if (iters <= 0) call abor1('ectrans_benchmark:iters <= 0')

!===================================================================================================

call ectrans_program_init(verbosity=verbosity, pinning=lpinning)

call ectrans_spectral_decomposition (nproc, nprtrw, nprtrv)
call ectrans_gridpoint_decomposition(nproc, nprgpns, nprgpew)

!===================================================================================================

call timer_init          %register_in_gstats('INIT           - Initialization')
call timer_setup_trans0  %register_in_gstats('SETUP_TRANS0   - Setup ecTrans')
call timer_setup_trans   %register_in_gstats('SETUP_TRANS    - Setup ecTrans handle')
call timer_loop          %register_in_gstats('LOOP           - Loop over iterations')
call timer_invtrans      %register_in_gstats('INV_TRANS      - Inverse transform')
call timer_dirtrans      %register_in_gstats('DIR_TRANS      - Direct transform')
call timer_tstep         %register_in_gstats('TIME STEP      - Time step')
call timer_norms         %register_in_gstats('NORMS          - Norm computation')
call timer_dump_checksums%register_in_gstats('DUMP_CHECKSUMS - Compute and dump checksums')
call timer_dump_values   %register_in_gstats('DUMP_VALUES    - Dump values for debugging')

!===================================================================================================

call timer_init%start()

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

if (verbosity >= 1) write(nout,'(a)')'======= Setup ecTrans ======='

call timer_setup_trans0%start()
call setup_trans0(kout=nout, kerr=nerr, kprintlev=merge(2, 0, verbosity == 1), kpromatr=npromatr, &
  &               kprgpns=nprgpns, kprgpew=nprgpew, kprtrw=nprtrw, ldsync_trans=lsync_trans,  &
  &               ldeq_regions=leq_regions, ldalloperm=lalloperm, ldmpoff=.not.ectrans_mpi_enabled(), &
  &               kopt_memory_tr=nopt_mem_tr)
call timer_setup_trans0%stop()

call timer_setup_trans%start()
call setup_trans(ksmax=nsmax, kdgl=ndgl, kloen=nloen, ldsplit=.true., lduserpnm=luserpnm, &
  &              lduseflt=luseflt)
call timer_setup_trans%stop()

if (ldump_checksums) then
  call invtrans_checksums_file_writer%open(trim(cchecksums_path)//'_inv_trans.checksums')
  call dirtrans_checksums_file_writer%open(trim(cchecksums_path)//'_dir_trans.checksums')
endif

!===================================================================================================
! Compute some dimensions and derived parameters
!===================================================================================================

call trans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg, kmysetv=mysetv)

if (nproma == 0) then ! no blocking (default when not specified)
  nproma = ngptot
endif

! Calculate number of NPROMA blocks
ngpblks = (ngptot - 1)/nproma+1

! Calculate the index of the last grid point in the last block
iend = ngptot - nproma * (ngpblks - 1)

! Determine the number of levels attributed to each member of the V set
numll = ectrans_make_spectral_distribution(nflevg, nprtrv)
nflevl = numll(mysetv)


!===================================================================================================
! Print information before starting
!===================================================================================================

! Print configuration details
if (verbosity >= 0 .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a)')'======= Start of runtime parameters ======='
  write(nout,'(" ")')
  write(nout,'("call_mode  ",i0)') icall_mode
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
  write(nout,'("npromatr   ",i0)') npromatr
  write(nout,'("ngpblks    ",i0)') ngpblks
  write(nout,'("nspec2     ",i0)') nspec2
  write(nout,'("nspec2g    ",i0)') nspec2g
  write(nout,'("luseflt    ",l1)') luseflt
  write(nout,'("nopt_mem_tr",i0)') nopt_mem_tr
  write(nout,'("lvordiv    ",l1)') lvordiv
  write(nout,'("lscders    ",l1)') lscders
  write(nout,'("luvder     ",l1)') luvder
  write(nout,'("lfield_api ",l1)') lfield_api
  write(nout,'("lalloperm  ",l1)') lalloperm
  write(nout,'(" ")')
  write(nout,'(a)') '======= End of runtime parameters ======='
  write(nout,'(" ")')
end if

!===================================================================================================
! Allocate and initialize spectral arrays
!===================================================================================================

! Compute spectral distribution variables for 3D fields
allocate(ivset(nflevg))
ilev = 0
do jb = 1, nprtrv
  do jlev=1, numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo

! Initialize vorticity and divergence - same for both call modes
call allocator%allocate('zspvor', zspvor, [nflevl,nspec2])
call allocator%allocate('zspdiv', zspdiv, [nflevl,nspec2])
call initialize_spectral_field(nsmax, zspvor)
call initialize_spectral_field(nsmax, zspdiv)

! Initialize spectral arrays differently depending on call mode
if (icall_mode == 1) then
  ! Compute spectral distribution variables for call mode 1's combined 2D/3D spectral array
  allocate(ivsetsc(nfld*nflevg+1))
  do ifld = 1, nfld
    ilev = 0
    do jb = 1, nprtrv
      do jlev = 1, numll(jb)
        ilev = ilev + 1
        ivsetsc(ilev + (ifld - 1)*nflevg) = jb
      enddo
    enddo
  enddo
  ivsetsc(nfld*nflevg+1) = 1

  call allocator%allocate('zspscalar', zspscalar, [count(ivsetsc == mysetv),nspec2])
  call initialize_spectral_field(nsmax, zspscalar)
else
  ivsetsc2(1) = min(nflevg+1, nprtrv)

  call allocator%allocate('zspsc3a', zspsc3a, [nflevl,nspec2,nfld])
  call allocator%allocate('zspsc2',  zspsc2,  [1,nspec2])
  do ifld = 1, nfld
    call initialize_spectral_field(nsmax, zspsc3a(:,:,ifld))
  enddo
  call initialize_spectral_field(nsmax, zspsc2)
endif

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
  call allocator%allocate('zgp', zgp, [nproma,itotal_fields,ngpblks])
else
  call allocator%allocate('zgpuv', zgpuv, [nproma,nflevg,inum_wind_fields,ngpblks])
  call allocator%allocate('zgp3a', zgp3a, [nproma,nflevg,inum_sc_3d_fields,ngpblks])
  call allocator%allocate('zgp2', zgp2, [nproma,inum_sc_2d_fields,ngpblks])
endif

!===================================================================================================
! Setup fields for field API if using
!===================================================================================================

#if USE_FIELD_API
if (lfield_api) then
  if (icall_mode == 1) then
    call wrap_benchmark_fields_zgp(ywflds, lvordiv, lscders, luvder, nflevg, 1 + nflevg * nfld, &
      &                            zspvor, zspdiv, zspscalar, zgp)
    call create_fields_lists(ywflds, ylf, kvsetuv=ivset, kvsetsc=ivsetsc)
  else
    call wrap_benchmark_fields(ywflds, lvordiv, lscders, luvder, 1, 1, nfld, zspvor, zspdiv, &
      &                        zspsc3a, zspsc2, zgpuv, zgp3a, zgp2)
    call create_fields_lists(ywflds, ylf, kvsetuv=ivset, kvsetsc2=ivsetsc2, kvsetsc=ivset)
  endif
endif
#endif

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
      write(nout,'("0x",Z16.16)') transfer(znormvor1(ifld),0_jpim)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv1(ifld)
      write(nout,'("0x",Z16.16)') transfer(znormdiv1(ifld),0_jpim)
    enddo
    if (icall_mode == 1) then
      do ifld = 1, nfld*nflevg+1
        write(nout,'("norm zspscalar(",i4,",:,1) = ",f20.15)') ifld, znormscalar1(ifld)
        write(nout,'("0x",Z16.16)') transfer(znormscalar1(ifld),0_jpim)
      enddo
    else
      if (nfld > 0) then
        do ifld = 1, nflevg
          write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormsc3a1(ifld)
          write(nout,'("0x",Z16.16)') transfer(znormsc3a1(ifld),0_jpim)
        enddo
      endif
      write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15)') 1, znormsc21(1)
      write(nout,'("0x",Z16.16)') transfer(znormsc21(1),0_jpim)
    endif
  endif
endif

if (verbosity >= 0 .and. myproc == 1) then
  write(nout,'(a,i0,a,f9.2,a)') "ectrans_benchmark initialisation, on ",nproc,&
                                & " tasks, took ",timer_init%elapsed()," sec"
endif

call timer_init%stop()

!===================================================================================================
! Do spectral transform loop
!===================================================================================================

if (verbosity >= 1 .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a)') '======= Start of spectral transforms  ======='
  write(nout,'(" ")')
endif

call ectrans_gstats_enable(.false.) ! Pause gstats timers during warmup iterations

write(nout,'(" ")')
write(nout,'(a,i0,a,i0,a)') 'Running for ', iters, ' iterations with ', iters_warmup, &
  & ' extra warm-up iterations (warmup marked with *)'
write(nout,'(" ")')

do jstep = 1, iters+iters_warmup
  if (jstep == iters_warmup + 1) then
    call ectrans_gstats_enable(.true.) ! Resume gstats timers after warmup iterations
    lloop_timer_active = .true.
    call timer_loop%start()
  endif

  call timer_tstep%start()

  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

  call timer_invtrans%start()
  if (lfield_api) then
#if USE_FIELD_API
    call inv_trans_field_api(kresol=1, ydfspscalar=ylf%spscalar, ydfspvor=ylf%spvor, &
      &                      ydfspdiv=ylf%spdiv, ydfscalar=ylf%scalar, ydfu=ylf%u, ydfv=ylf%v, &
      &                      ydfvor=ylf%vor, ydfdiv=ylf%div, ydfscalar_ns=ylf%scalar_ns, &
      &                      ydfscalar_ew=ylf%scalar_ew, ydfu_ew=ylf%u_ew, ydfv_ew=ylf%v_ew, &
      &                      kgptot = ngptot)
    call synchost_rdonly_wrapped_fields(ywflds)
#else
    call abor1('ectrans_benchmark: No field API support')
#endif
  else if (icall_mode == 1) then
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
  call timer_invtrans%stop()

  if (ldump_checksums .and. (.not. lloop_timer_active .or. iters_warmup == 0)) then
    call ectrans_gstats_enable(.true.)
    call timer_dump_checksums%start()
    call invtrans_checksums_file_writer%write_iteration_separator(jstep)
    if (icall_mode == 1) then
      ! Remove trash at end of last block for checksumming
      zgp (iend+1:, :, ngpblks) = 0
      call invtrans_checksums_file_writer%write_checksums_pgp(zgp)
    else
      ! Remove trash at end of last block for checksumming
      zgpuv (iend+1:, :, :, ngpblks) = 0
      zgp3a (iend+1:, :, :, ngpblks) = 0
      zgp2  (iend+1:,    :, ngpblks) = 0
      call invtrans_checksums_file_writer%write_checksums_pgp_uv_3a_2(zgpuv=zgpuv, zgp3a=zgp3a, zgp2=zgp2)
    endif
    call timer_dump_checksums%stop()
    call ectrans_gstats_enable(lloop_timer_active)
  endif

  if (lloop_timer_active) then
    call timings_invtrans%push_back(timer_invtrans%elapsed())
  endif

  !=================================================================================================
  ! While in grid point space, dump the values to disk, for debugging only
  !=================================================================================================

  if (ldump_values .and. mod(jstep,10) == 1) then
    call timer_dump_values%start()
    ! dump a field to a binary file
    if (myproc == 1) then
      allocate(global_field(ngptotg,1))
    endif
    if (icall_mode == 1) then
      islice = (ipgpuv_end - 1) * nflevg
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp(:,islice:islice,:), 'U')
      islice = ipgpuv_end * nflevg
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp(:,islice:islice,:), 'V')
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp(:,ipgp_end:ipgp_end,:), 'S')
      islice = ipgp_end - 1
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp(:,islice:islice,:), 'T')
    else
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgpuv(:,nflevg:nflevg,1,:), 'U')
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgpuv(:,nflevg:nflevg,2,:), 'V')
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp2(:,1:1,:), 'S')
      call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp3a(:,nflevg:nflevg,1,:), 'T')
    endif
    if (myproc == 1) then
      deallocate(global_field)
    endif
    call timer_dump_values%stop()
  endif

  !=================================================================================================
  ! Do direct transform
  !=================================================================================================

  call timer_dirtrans%start()
  if (lfield_api) then
#if USE_FIELD_API
    call dir_trans_field_api(kresol=1, ydfscalar=ylf%scalar, ydfu=ylf%u, ydfv=ylf%v, &
      &                      ydfspscalar=ylf%spscalar, ydfspvor=ylf%spvor, ydfspdiv=ylf%spdiv)
    call synchost_rdonly_wrapped_fields(ywflds)
#else
    call abor1('ectrans_benchmark: No field API support')
#endif
  else if (icall_mode == 1) then
    call dir_trans(pgp=zgp(:,ipgp_start:ipgp_end,:), pspvor=zspvor, pspdiv=zspdiv, &
      &            pspscalar=zspscalar, kvsetuv=ivset, kvsetsc=ivsetsc, kproma=nproma)
  else
    call dir_trans(pgpuv=zgpuv(:,:,ipgpuv_start:ipgpuv_end,:), &
      &            pgp3a=zgp3a(:,:,1:nfld,:), pgp2=zgp2(:,1:1,:), &
      &            pspvor=zspvor, pspdiv=zspdiv, pspsc3a=zspsc3a, pspsc2=zspsc2, &
      &            kvsetuv=ivset, kvsetsc2=ivsetsc2, kvsetsc3a=ivset, kproma=nproma)
  endif
  call timer_dirtrans%stop()

  if (ldump_checksums .and. (.not. lloop_timer_active .or. iters_warmup == 0)) then
    call ectrans_gstats_enable(.true.)
    call timer_dump_checksums%start()
    call dirtrans_checksums_file_writer%write_iteration_separator(jstep)
    if (icall_mode == 1) then
      call dirtrans_checksums_file_writer%write_checksums_psp(ivset=ivset, ivsetsc=ivsetsc, &
        &                     zspvor=zspvor, zspdiv=zspdiv, zspscalar=zspscalar)
    else
      call dirtrans_checksums_file_writer%write_checksums_psp_3a_2(ivset=ivset, ivsetsc2=ivsetsc2, &
        &                          zspvor=zspvor, zspdiv=zspdiv, zspsc3a=zspsc3a, zspsc2=zspsc2)
    endif
    call timer_dump_checksums%stop()
    call ectrans_gstats_enable(lloop_timer_active)
  endif
  call timer_tstep%stop()
  if (lloop_timer_active) then
    call timings_dirtrans%push_back(timer_dirtrans%elapsed())
    call timings_tstep%push_back(timer_tstep%elapsed())
  endif

  !=================================================================================================
  ! Print norms
  !=================================================================================================

  if (lprint_norms) then
    call timer_norms%start()
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
        if (lloop_timer_active) then
          write(nout,'("time step ",i6," took", f8.4," | zspvor max err=",e10.3,&
          & " | zspdiv max err=",e10.3," | zspscalar max err=",e10.3)') &
          &  jstep, timer_tstep%elapsed(), zmaxerr(1), zmaxerr(2), zmaxerr(3)
        else
          write(nout,'("time step *",i5," took", f8.4," | zspvor max err=",e10.3,&
          & " | zspdiv max err=",e10.3," | zspscalar max err=",e10.3)') &
          &  jstep, timer_tstep%elapsed(), zmaxerr(1), zmaxerr(2), zmaxerr(3)
        endif
      else
        zmaxerr(4) = maxval(abs((znormsc21 / znormsc2) - 1.0_jprb))
        if (nfld > 0) then
          zmaxerr(3) = maxval(abs((znormsc3a1 / znormsc3a) - 1.0_jprb))
          if (lloop_timer_active) then
            write(nout,'("time step ",i6," took", f8.4," | zspvor max err=",e10.3,&
            & " | zspdiv max err=",e10.3," | zspsc3a max err=",e10.3," | zspsc2 max err=",e10.3)') &
            &  jstep, timer_tstep%elapsed(), zmaxerr(1), zmaxerr(2), zmaxerr(3), zmaxerr(4)
          else
            write(nout,'("time step *",i5," took", f8.4," | zspvor max err=",e10.3,&
            & " | zspdiv max err=",e10.3," | zspsc3a max err=",e10.3," | zspsc2 max err=",e10.3)') &
            &  jstep, timer_tstep%elapsed(), zmaxerr(1), zmaxerr(2), zmaxerr(3), zmaxerr(4)
          endif
        else
          if (lloop_timer_active) then
            write(nout,'("time step ",i6," took", f8.4," | zspvor max err=",e10.3,&
                        & " | zspdiv max err=",e10.3," | zspsc2 max err=",e10.3)') &
                        &  jstep, timer_tstep%elapsed(), zmaxerr(1), zmaxerr(2), zmaxerr(4)
          else
            write(nout,'("time step *",i5," took", f8.4," | zspvor max err=",e10.3,&
                        & " | zspdiv max err=",e10.3," | zspsc2 max err=",e10.3)') &
                        &  jstep, timer_tstep%elapsed(), zmaxerr(1), zmaxerr(2), zmaxerr(4)
          endif
        endif
      endif
    endif
    call timer_norms%stop()
  else
    if (lloop_timer_active) then
      write(nout,'("time step ",i6," took", f8.4)') jstep, timer_tstep%elapsed()
    else
      write(nout,'("time step *",i5," took", f8.4)') jstep, timer_tstep%elapsed()
    endif
  endif
enddo

!===================================================================================================

call timer_loop%stop()
call timings_loop%push_back(timer_loop%elapsed())

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
        write(nout,'("0x",Z16.16)') transfer(znormvor(ifld), 0_jpim)
      enddo
    endif
    zmaxerr(2) = maxval(abs((real(znormdiv1,jprd) / (real(znormdiv,jprd)) - 1.0_jprd)))
    if (verbosity >= 1) then
      do ifld = 1, nflevg
        write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv(ifld)
        write(nout,'("0x",Z16.16)') transfer(znormdiv(ifld), 0_jpim)
      enddo
    endif
    if (icall_mode == 1) then
      zmaxerr(3) = maxval(abs((znormscalar1 / znormscalar) - 1.0_jprb))
      if (verbosity >= 1) then
        do ifld = 1, nfld*nflevg+1
          write(nout,'("norm znormscalar( ",i4,",:)   = ",f20.15)') ifld, znormscalar(ifld)
          write(nout,'("0x",Z16.16)') transfer(znormscalar(ifld), 0_jpim)
        enddo
      endif
    else
      zmaxerr(4) = maxval(abs((znormsc21 / znormsc2) - 1.0_jprb))
      if (verbosity >= 1) then
        write(nout,'("norm znormsc2( ",i4,",:)   = ",f20.15)') 1, znormsc2(1)
        write(nout,'("0x",Z16.16)') transfer(znormsc2(1), 0_jpim)
      endif
      if (nfld > 0) then
        zmaxerr(3) = maxval(abs((znormsc3a1 / znormsc3a) - 1.0_jprb))
        if (verbosity >= 1) then
          do ifld = 1, nflevg
            write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormsc3a(ifld)
            write(nout,'("0x",Z16.16)') transfer(znormsc3a(ifld), 0_jpim)
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
    if (ectrans_mpi_enabled()) then
      call mpl_broadcast(ierr,kroot=1,ktag=1)
    endif

    ! Halt if correctness checker failed
    if (ierr == 1) then
      error stop
    endif
  endif
endif

!===================================================================================================
! Report timings
!===================================================================================================

write(nout,'(a)') '======= Start of time step stats ======='
write(nout,'(" ")')
write(nout,'("Inverse transforms")')
write(nout,'("------------------")')
write(nout,'("avg  (s): ",f8.4)') timings_invtrans%global_avg()
write(nout,'("min  (s): ",f8.4)') timings_invtrans%global_min()
write(nout,'("max  (s): ",f8.4)') timings_invtrans%global_max()
write(nout,'("med  (s): ",f8.4)') timings_invtrans%global_med()
write(nout,'(" ")')
write(nout,'("Direct transforms")')
write(nout,'("-----------------")')
write(nout,'("avg  (s): ",f8.4)') timings_dirtrans%global_avg()
write(nout,'("min  (s): ",f8.4)') timings_dirtrans%global_min()
write(nout,'("max  (s): ",f8.4)') timings_dirtrans%global_max()
write(nout,'("med  (s): ",f8.4)') timings_dirtrans%global_med()
write(nout,'(" ")')
write(nout,'("Inverse-direct transforms")')
write(nout,'("-------------------------")')
write(nout,'("avg  (s): ",f8.4)') timings_tstep%global_avg()
write(nout,'("min  (s): ",f8.4)') timings_tstep%global_min()
write(nout,'("max  (s): ",f8.4)') timings_tstep%global_max()
write(nout,'("med  (s): ",f8.4)') timings_tstep%global_med()
write(nout,'(" ")')
write(nout,'("loop (s): ",f8.4)') timings_loop%global_avg()
write(nout,'(" ")')
write(nout,'(a)') '======= End of time step stats ======='

!===================================================================================================
! Cleanup
!===================================================================================================
#if USE_FIELD_API
if (lfield_api) then
  call delete_wrapped_fields(ywflds)
  call delete_fields_lists(ylf)
endif
#endif

call allocator%deallocate('zspvor', zspvor)
call allocator%deallocate('zspdiv', zspdiv)

if (icall_mode == 1) then
  call allocator%deallocate('zspscalar', zspscalar)
else
  call allocator%deallocate('zspsc3a', zspsc3a)
  call allocator%deallocate('zspsc2', zspsc2)
endif

if (icall_mode == 1) then
  call allocator%deallocate('zgp', zgp)
else
  call allocator%deallocate('zgpuv', zgpuv)
  call allocator%deallocate('zgp3a', zgp3a)
  call allocator%deallocate('zgp2', zgp2)
endif

!===================================================================================================

if (lmeminfo) then
  write(nout,*)
  call ectrans_print_memory_usage()
endif

if (ldump_checksums) then
  if (verbosity >= 1) write(nout,*)
  call invtrans_checksums_file_writer%close()
  call dirtrans_checksums_file_writer%close()
endif

write(nout,*)
write(nout,'(A)') 'ectrans benchmark completed successfully'
call ectrans_program_end()

!===================================================================================================

contains

!===================================================================================================

subroutine print_help(unit)
  integer(kind=jpim), intent(in) :: unit
  integer :: nout
  nout = unit

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
  write(nout, "(a)") "    --npromatr NPROMATR Perform transforms in blocks of size NPROMATR rather&
    & than all at once"
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
  write(nout, "(a)") "    --no-pinning        Disable memory-pinning (a.k.a. page-locked memory) &
   & to allocate fields for GPU version"
  write(nout, "(a)") "    --field-api         Use the field api interface of ecTrans"
  write(nout, "(a)") "    --callmode          The call mode for INV_TRANS and DIR_TRANS (1 or 2)"
  write(nout, "(a)") "                        Call mode 1 uses arrays PSPVOR, PSPDIV, PSPSCALAR and&
   & PGP"
  write(nout, "(a)") "                        Call mode 2 uses arrays PSPVOR, PSPDIV, PSPSC3A,&
   & PSPSC3B, PSPSC2, PGPUV, PGP3A, PGP3B, PGP2"
  write(nout, "(a)") "                        See&
   & https://sites.ecmwf.int/docs/ectrans/page/api.html for more information (default  = 2)"
  write(nout, "(a)") "    --deallocate-foubuf-temps Enable deallocation of temporary Fourier-space&
   & buffers (default = off, when enabled equivalent to LALLOPERM=.FALSE.)"
  write(nout, "(a)") ""
  write(nout, "(a)") "DEBUGGING"
  write(nout, "(a)") "    --dump-values             Output gridpoint fields in unformatted binary file"
  write(nout, "(a)") "    --dump-checksums FILENAME Output CRC64 checksums of fields in text file named FILENAME"
  write(nout, "(a)") "                              Checksums are written only during warmup iterations (see --niter-warmup)"
  write(nout, "(a)") ""

end subroutine print_help

!===================================================================================================

subroutine get_command_line_arguments(nsmax, cgrid, iters, iters_warmup, nfld, nlev, lvordiv, &
  &                                   lscders, luvder, luseflt, nopt_mem_tr, nproma, npromatr, &
  &                                   verbosity, ldump_values, lprint_norms, lmeminfo, nprtrv, &
  &                                   nprtrw, ncheck, lpinning, lfield_api, icall_mode, ldump_checksums, &
  &                                   cchecksums_path, lalloperm)

  integer, intent(inout) :: nsmax           ! Spectral truncation
  character(len=16), intent(inout) :: cgrid ! Grid
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
  integer, intent(inout) :: npromatr        ! block size for field-blocking
  integer, intent(inout) :: verbosity       ! Level of verbosity
  logical, intent(inout) :: ldump_values    ! Dump values of grid point fields for debugging
  logical, intent(inout) :: ldump_checksums ! Dump CRC checksums
  logical, intent(inout) :: lprint_norms    ! Calculate and print spectral norms of fields
  logical, intent(inout) :: lmeminfo        ! Show information from FIAT ec_meminfo routine at the
                                            ! end
  integer, intent(inout) :: nprtrv          ! Size of V set (spectral decomposition)
  integer, intent(inout) :: nprtrw          ! Size of W set (spectral decomposition)
  integer, intent(inout) :: ncheck          ! The multiplier of the machine epsilon used as a
                                            ! tolerance for correctness checking
  logical, intent(inout) :: lpinning        ! Use memory-pinning (a.k.a. page-locked memory) to allocate fields for GPU version
  logical, intent(inout) :: lfield_api
  integer, intent(inout) :: icall_mode      ! The call mode for inv_trans and dir_trans
                                            ! 1: pspvor, pspdiv, pspscalar, pgp
                                            ! 2: pspvor, pspdiv, pspsc3a, pspsc2, pgpuv, pgp3a, pgp2

  character(len=128), intent(inout) :: cchecksums_path ! path to export checksum files
  logical, intent(inout) :: lalloperm                  ! keep FOUBUF & FOUBUF_IN allocated
  character(len=128) :: carg          ! Storage variable for command line arguments
  type(ectrans_command_line_parser) :: parser

  call parser%setup(print_help)
  do while (parser%next_arg(carg))
    select case(carg)
      ! Parse help argument
      case('-h', '--help')
        call parser%print_help()
        stop
      ! Parse verbosity argument
      case('-v')
        verbosity = 1
      ! Parse number of iterations argument
      case('-n', '--niter')
        iters = parser%get_int_value()
        if (iters < 1) then
          call parser%parsing_failed("Invalid argument for -n: must be > 0")
        end if
      case('--niter-warmup')
        iters_warmup = parser%get_int_value()
        if (iters_warmup < 0) then
          call parser%parsing_failed("Invalid argument for --niter-warmup: must be >= 0")
        end if
      ! Parse spectral truncation argument
      case('-t', '--truncation')
        nsmax = parser%get_int_value()
        if (nsmax < 1) then
          call parser%parsing_failed("Invalid argument for -t: must be > 0")
        end if
      case('-g', '--grid'); cgrid = parser%get_str_value()
      case('-f', '--nfld'); nfld = parser%get_int_value()
      case('-l', '--nlev'); nlev = parser%get_int_value()
      case('--vordiv'); lvordiv = .true.
      case('--scders'); lscders = .true.
      case('--uvders'); luvder = .true.
      case('--flt'); luseflt = .true.
      case('--mem-tr'); nopt_mem_tr = parser%get_int_value()
      case('--nproma'); nproma = parser%get_int_value()
      case('--npromatr'); npromatr = parser%get_int_value()
      case('--dump-values'); ldump_values = .true.
      case('--dump-checksums')
        ldump_checksums = .true.
        cchecksums_path = parser%get_str_value()
      case('--norms'); lprint_norms = .true.
      case('--meminfo'); lmeminfo = .true.
      case('--nprtrv'); nprtrv = parser%get_int_value()
      case('--nprtrw'); nprtrw = parser%get_int_value()
      case('-c', '--check'); ncheck = parser%get_int_value()
      case('--no-pinning'); lpinning = .False.
      case('--field-api'); lfield_api = .True.
      case('--callmode')
          icall_mode = parser%get_int_value()
          if (icall_mode < 1 .or. icall_mode > 2) then
            call parser%parsing_failed("Invalid argument for --callmode: must be 1 or 2")
          end if
      case('--deallocate-foubuf-temps'); lalloperm = .false.
      case default
        call parser%parsing_failed("Unrecognised argument: " // trim(carg))
    end select
  end do

end subroutine get_command_line_arguments

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

  integer :: num_my_zon_wns
  integer, allocatable :: my_zon_wns(:)

  ! Choose a spherical harmonic to initialize arrays
  integer, parameter :: m_num = 4  ! Zonal wavenumber
  integer, parameter :: l_num = 19  ! Total wavenumber

  ! First initialise all spectral coefficients to zero
  field(:) = 0.0

  ! Get zonal wavenumbers this rank is responsible for
  call trans_inq(knump=num_my_zon_wns)
  allocate(my_zon_wns(num_my_zon_wns))
  call trans_inq(kmyms=my_zon_wns)

  ! If rank is responsible for the chosen zonal wavenumber...
  if (any(my_zon_wns == m_num) ) then
    block
      integer, allocatable :: nasm0(:)
      integer :: index

      ! Get array of spectral array addresses (this maps (m, n=m) to array index)
      allocate(nasm0(0:nsmax))
      call trans_inq(kasm0=nasm0)

      ! Find out local array index of chosen spherical harmonic
      index = nasm0(m_num) + 2 * (l_num - m_num) + 1

      ! Set just that element to a constant value
      field(index) = 1.0
    end block
  end if

end subroutine initialize_2d_spectral_field

!===================================================================================================

subroutine dump_gridpoint_field(jstep, myproc, nproma, gfld, fld, fldchar)
  use parkind1, only : jprb, jpim
  implicit none
  ! Dump a 2d field to a binary file.

  integer(kind=jpim), intent(in) :: jstep ! Time step, used for naming file
  integer(kind=jpim), intent(in) :: myproc ! MPI rank, used for naming file
  integer(kind=jpim), intent(in) :: nproma ! Size of nproma
  real(kind=jprb)   , intent(inout) :: gfld(:,:) ! 2d global field
  real(kind=jprb)   , intent(in) :: fld(:,:,:) ! 3d local field
  character         , intent(in) :: fldchar ! Single character field identifier

  character(len=10)  :: filename
  integer(kind=jpim) :: fileunit ! unit number for output file
  integer(kind=jpim) :: ilev

#include "gath_grid.h"

  filename = "x.xxxx.dat"
  if (myproc == 1) then
    write(filename(1:1),'(a1)') fldchar
    write(filename(3:6),'(i4.4)') jstep
    open(newunit=fileunit, file=filename, form='unformatted')
  endif
  do ilev=1,size(fld,2)
    call gath_grid(gfld(:,:),nproma,1,(/1/),1,fld(:,ilev:ilev,:))
    if (myproc == 1) write(unit=fileunit) gfld(:,1)
  enddo
  if (myproc == 1) then
    close(fileunit)
  endif

end subroutine dump_gridpoint_field

!===================================================================================================

end program ectrans_benchmark

!===================================================================================================
