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
!

use parkind1, only: jpim, jprd
use oml_mod ,only : oml_max_threads
use mpl_module
use yomgstats, only: jpmaxstat
use yomhook, only : dr_hook_init
use transform_driver_mod_dp, only :& 
 ectrans_setup_dp =>     ectrans_setup, & 
 ectrans_setup0_dp =>     ectrans_setup0, &
 ectrans_trans_inq_dp =>     ectrans_trans_inq, &
 ectrans_allocate_spectral_dp =>     ectrans_allocate_spectral, &
 ectrans_allocate_grid_dp =>     ectrans_allocate_grid, &
 ectrans_deallocate_grid_dp =>     ectrans_deallocate_grid, &
 ectrans_allocate_normdata_dp =>     ectrans_allocate_normdata,  &
 ectrans_calculate_norms_dp =>     ectrans_calculate_norms, & 
 ectrans_print_norms_init_dp =>     ectrans_print_norms_init, &
 ectrans_allocate_timers_dp =>     ectrans_allocate_timers, &
 ectrans_set_ztstep_start_dp =>     ectrans_set_ztstep_start, &
 ectrans_set_ztstep_end_dp =>     ectrans_set_ztstep_end, & 
 ectrans_inv_trans_dp =>     ectrans_inv_trans,&
 ectrans_dump_dp =>     ectrans_dump, &
 ectrans_direct_trans_dp =>     ectrans_direct_trans, &
 ectrans_calculate_timings_dp =>     ectrans_calculate_timings, &
 ectrans_print_norms_calc_dp =>     ectrans_print_norms_calc, &
 ectrans_print_norms_fin_dp =>     ectrans_print_norms_fin, &
 ectrans_print_norms_fails_dp =>     ectrans_print_norms_fails, &
 ectrans_norms_reduce_dp =>     ectrans_norms_reduce, &
 ectrans_compute_time_stats_dp =>     ectrans_compute_time_stats, &
 ectrans_print_time_stats_dp =>     ectrans_print_time_stats, &
 ectrans_print_timestep_dp =>     ectrans_print_timestep

use transform_driver_mod_sp, only :& 
 ectrans_setup_sp =>     ectrans_setup, & 
 ectrans_setup0_sp =>     ectrans_setup0, &
 ectrans_trans_inq_sp =>     ectrans_trans_inq, &
 ectrans_allocate_spectral_sp =>     ectrans_allocate_spectral, &
 ectrans_allocate_grid_sp =>     ectrans_allocate_grid, &
 ectrans_deallocate_grid_sp =>     ectrans_deallocate_grid, &
 ectrans_allocate_normdata_sp =>     ectrans_allocate_normdata,  &
 ectrans_calculate_norms_sp =>     ectrans_calculate_norms, & 
 ectrans_print_norms_init_sp =>     ectrans_print_norms_init, &
 ectrans_allocate_timers_sp =>     ectrans_allocate_timers, &
 ectrans_set_ztstep_start_sp =>     ectrans_set_ztstep_start, &
 ectrans_set_ztstep_end_sp =>     ectrans_set_ztstep_end, & 
 ectrans_inv_trans_sp =>     ectrans_inv_trans,&
 ectrans_dump_sp =>     ectrans_dump, &
 ectrans_direct_trans_sp =>     ectrans_direct_trans, &
 ectrans_calculate_timings_sp =>     ectrans_calculate_timings, &
 ectrans_print_norms_calc_sp =>     ectrans_print_norms_calc, &
 ectrans_print_norms_fin_sp =>     ectrans_print_norms_fin, &
 ectrans_print_norms_fails_sp =>     ectrans_print_norms_fails, &
 ectrans_norms_reduce_sp =>     ectrans_norms_reduce, &
 ectrans_compute_time_stats_sp =>     ectrans_compute_time_stats, &
 ectrans_print_time_stats_sp =>     ectrans_print_time_stats, &
 ectrans_print_timestep_sp =>     ectrans_print_timestep
implicit none
 procedure(ectrans_setup_dp),pointer ::     ectrans_setup  
 procedure(ectrans_setup0_dp),pointer ::     ectrans_setup0 
 procedure(ectrans_trans_inq_dp),pointer ::     ectrans_trans_inq 
 procedure(ectrans_allocate_spectral_dp ),pointer ::     ectrans_allocate_spectral 
 procedure(ectrans_allocate_grid_dp),pointer ::     ectrans_allocate_grid 
 procedure(ectrans_deallocate_grid_dp),pointer ::     ectrans_deallocate_grid 
 procedure(ectrans_allocate_normdata_dp),pointer ::     ectrans_allocate_normdata
 procedure(ectrans_calculate_norms_dp),pointer ::     ectrans_calculate_norms  
 procedure(ectrans_print_norms_init_dp),pointer ::     ectrans_print_norms_init 
 procedure(ectrans_allocate_timers_dp),pointer ::     ectrans_allocate_timers 
 procedure(ectrans_set_ztstep_start_dp),pointer ::     ectrans_set_ztstep_start 
 procedure(ectrans_set_ztstep_end_dp),pointer ::     ectrans_set_ztstep_end  
 procedure(ectrans_inv_trans_dp),pointer ::     ectrans_inv_trans
 procedure(ectrans_dump_dp),pointer ::     ectrans_dump 
 procedure(ectrans_direct_trans_dp),pointer ::     ectrans_direct_trans 
 procedure(ectrans_calculate_timings_dp),pointer ::     ectrans_calculate_timings 
 procedure(ectrans_print_norms_calc_dp),pointer ::     ectrans_print_norms_calc 
 procedure(ectrans_print_norms_fin_dp),pointer ::     ectrans_print_norms_fin 
 procedure(ectrans_print_norms_fails_dp),pointer ::     ectrans_print_norms_fails 
 procedure(ectrans_norms_reduce_dp),pointer ::     ectrans_norms_reduce 
 procedure(ectrans_compute_time_stats_dp),pointer ::     ectrans_compute_time_stats 
 procedure(ectrans_print_time_stats_dp),pointer ::     ectrans_print_time_stats 
 procedure(ectrans_print_timestep_dp),pointer ::     ectrans_print_timestep

integer(kind=jpim) :: istack, getstackusage
!
! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT
integer(kind=jpim), parameter :: noutdump = 7 ! Unit number for field output
!
!! Default parameters
integer(kind=jpim) :: nsmax   = 79  ! Spectral truncation
integer(kind=jpim) :: iters   = 10  ! Number of iterations for transform test
integer(kind=jpim) :: nfld    = 1   ! Number of scalar fields 
integer(kind=jpim) :: nlev    = 1   ! Number of vertical levels
!
integer(kind=jpim) :: nflevg
integer(kind=jpim) :: ndgl ! Number of latitudes
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
!
integer(kind=jpim), allocatable :: nloen(:), nprcids(:)
integer(kind=jpim) :: myproc, jj
integer :: jstep
!
real(kind=jprd) :: ztinit, ztloop, timef 
real(kind=jprd) :: zaveave(0:jpmaxstat)
logical :: lstack = .false. ! Output stack info
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
logical :: lvordiv = .false.
logical :: lscders = .false.
logical :: luvders = .false.
logical :: lprint_norms = .true. ! Calculate and print spectral norms
logical :: lmeminfo = .false. ! Show information from FIAT routine ec_meminfo at the end
!
integer(kind=jpim) :: nstats_mem = 0
integer(kind=jpim) :: ntrace_stats = 0
integer(kind=jpim) :: nprnt_stats = 1
!
! The multiplier of the machine epsilon used as a tolerance for correctness checking
! ncheck = 0 (the default) means that correctness checking is disabled
integer(kind=jpim) :: ncheck = 0
!
logical :: lmpoff = .false. ! Message passing switch
!
! Verbosity level (0 or 1)
integer :: verbosity = 0
!
!integer(kind=jpim) :: ncombflen = 1800000 ! Size of comm buffer
!
integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nthread
integer(kind=jpim) :: nprgpns ! Grid-point decomp
integer(kind=jpim) :: nprgpew ! Grid-point decomp
integer(kind=jpim) :: nprtrv = 0 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
integer(kind=jpim) :: nspecresmin = 80 ! Minimum spectral resolution, for controlling nprtrw
integer(kind=jpim) :: mysetv
integer(kind=jpim) :: mysetw
integer(kind=jpim) :: mp_type = 2 ! Message passing type
integer(kind=jpim) :: mbx_size = 150000000 ! Mailbox size
!
integer(kind=jpim), allocatable :: numll(:), ivset(:)
integer(kind=jpim) :: ivsetsc(1)
!
integer(kind=jpim) :: nflevl
!
! sumpini
integer(kind=jpim) :: isqr
!
integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: ngpblks
!! locals
integer(kind=jpim) :: iprtrv
integer(kind=jpim) :: iprtrw
integer(kind=jpim) :: iprused, ilevpp, irest, ilev, jlev
!
logical :: ldump_values = .false.
!
integer, external :: ec_mpirank
logical :: luse_mpi = .true.
!
character(len=16) :: cgrid = ''
!
integer(kind=jpim) :: ierr
integer(kind=jpim) :: ibackend

!===================================================================================================

#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"
#include "ec_meminfo.intfb.h"

!===================================================================================================

luse_mpi = detect_mpirun()

! Setup
call get_command_line_arguments(nsmax, cgrid, iters, nfld, nlev, lvordiv, lscders, luvders, &
  & luseflt, nproma, verbosity, ldump_values, lprint_norms, lmeminfo, nprtrv, nprtrw, ncheck)
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
  lstats_mpl    = .true.
  lstatscpu     = .true.
  nprnt_stats   = nproc
!  lstats_mem   = .true.
!  lstats_alloc = .true.
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
isqr = int(sqrt(real(nproc,jprd)))
do ja = isqr, nproc
  ib = nproc/ja
  if (ja*ib == nproc) then
    nprgpns = max(ja,ib)
    nprgpew = min(ja,ib)
    exit
  endif
enddo

! From sumpini, although this should be specified in namelist
if (nspecresmin == 0) nspecresmin = nproc

! Compute nprtrv and nprtrw if not provided on the command line
if (nprtrv > 0 .or. nprtrw > 0) then
  if (nprtrv == 0) nprtrv = nproc/nprtrw
  if (nprtrw == 0) nprtrw = nproc/nprtrv
  if (nprtrw*nprtrv /= nproc) call abor1('transform_test:nprtrw*nprtrv /= nproc')
  if (nprtrw > nspecresmin) call abor1('transform_test:nprtrw > nspecresmin')
else
  do jprtrv = 4, nproc
    nprtrv = jprtrv
    nprtrw = nproc/nprtrv
    if (nprtrv*nprtrw /= nproc) cycle
    if (nprtrv > nprtrw) exit
    if (nprtrw > nspecresmin) cycle
    if (nprtrw <= nspecresmin/(2*oml_max_threads())) exit
  enddo
  ! Go for approx square partition for backup
  if (nprtrv*nprtrw /= nproc .or. nprtrw > nspecresmin .or. nprtrv > nprtrw) then
    isqr = int(sqrt(real(nproc,jprd)))
    do ja = isqr, nproc
      ib = nproc/ja
      if (ja*ib == nproc) then
        nprtrw = max(ja, ib)
        nprtrv = min(ja, ib)
        if (nprtrw > nspecresmin ) then
          call abor1('transform_test:nprtrw (approx square value) > nspecresmin')
        endif
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
    call abor1('transform_test:inconsistency when computing mysetw and mysetv')
  endif
endif

if (.not. lmpoff) then
  call mpl_buffer_method(kmp_type=mp_type, kmbx_size=mbx_size, kprocids=nprcids, ldinfo=(verbosity>=1))
endif

! Determine number of local levels for fourier and legendre calculations
! based on the values of nflevg and nprtrv
allocate(numll(nprtrv+1))

! Calculate remainder
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
numll(iprused+1:nprtrv+1) = 0

nflevl = numll(mysetv)

ivsetsc(1) = iprused

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
  call gstats_labels(0)
endif
do ibackend=1,2
if(ibackend == 1) write(nout,'(a)')'======= Running benchmark #1 (sp) ======='
if(ibackend == 2) write(nout,'(a)')'======= Running benchmark #2 (dp) ======='
if(ibackend == 1) then
 ectrans_setup =>     ectrans_setup_sp  
 ectrans_setup0 =>     ectrans_setup0_sp 
 ectrans_trans_inq =>     ectrans_trans_inq_sp 
 ectrans_allocate_spectral =>     ectrans_allocate_spectral_sp 
 ectrans_allocate_grid =>     ectrans_allocate_grid_sp 
 ectrans_deallocate_grid =>     ectrans_deallocate_grid_sp 
 ectrans_allocate_normdata =>     ectrans_allocate_normdata_sp  
 ectrans_calculate_norms =>     ectrans_calculate_norms_sp  
 ectrans_print_norms_init =>     ectrans_print_norms_init_sp 
 ectrans_allocate_timers =>     ectrans_allocate_timers_sp 
 ectrans_set_ztstep_start =>     ectrans_set_ztstep_start_sp 
 ectrans_set_ztstep_end =>     ectrans_set_ztstep_end_sp  
 ectrans_inv_trans =>     ectrans_inv_trans_sp 
 ectrans_dump =>     ectrans_dump_sp 
 ectrans_direct_trans =>     ectrans_direct_trans_sp 
 ectrans_calculate_timings =>     ectrans_calculate_timings_sp 
 ectrans_print_norms_calc =>     ectrans_print_norms_calc_sp 
 ectrans_print_norms_fin =>     ectrans_print_norms_fin_sp 
 ectrans_print_norms_fails =>     ectrans_print_norms_fails_sp 
 ectrans_norms_reduce =>     ectrans_norms_reduce_sp 
 ectrans_compute_time_stats =>     ectrans_compute_time_stats_sp 
 ectrans_print_time_stats =>     ectrans_print_time_stats_sp 
 ectrans_print_timestep =>     ectrans_print_timestep_sp

elseif(ibackend == 2) then
 ectrans_setup =>     ectrans_setup_dp  
 ectrans_setup0 =>     ectrans_setup0_dp 
 ectrans_trans_inq =>     ectrans_trans_inq_dp 
 ectrans_allocate_spectral =>     ectrans_allocate_spectral_dp 
 ectrans_allocate_grid =>     ectrans_allocate_grid_dp 
 ectrans_deallocate_grid =>     ectrans_deallocate_grid_dp 
 ectrans_allocate_normdata =>     ectrans_allocate_normdata_dp  
 ectrans_calculate_norms =>     ectrans_calculate_norms_dp  
 ectrans_print_norms_init =>     ectrans_print_norms_init_dp 
 ectrans_allocate_timers =>     ectrans_allocate_timers_dp 
 ectrans_set_ztstep_start =>     ectrans_set_ztstep_start_dp 
 ectrans_set_ztstep_end =>     ectrans_set_ztstep_end_dp  
 ectrans_inv_trans =>     ectrans_inv_trans_dp 
 ectrans_dump =>     ectrans_dump_dp 
 ectrans_direct_trans =>     ectrans_direct_trans_dp 
 ectrans_calculate_timings =>     ectrans_calculate_timings_dp 
 ectrans_print_norms_calc =>     ectrans_print_norms_calc_dp 
 ectrans_print_norms_fin =>     ectrans_print_norms_fin_dp 
 ectrans_print_norms_fails =>     ectrans_print_norms_fails_dp 
 ectrans_norms_reduce =>     ectrans_norms_reduce_dp 
 ectrans_compute_time_stats =>     ectrans_compute_time_stats_dp 
 ectrans_print_time_stats =>     ectrans_print_time_stats_dp 
 ectrans_print_timestep =>     ectrans_print_timestep_dp
endif !backend 
!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

if (verbosity >= 1) write(nout,'(a)')'======= Setup ecTrans ======='

call gstats(1, 0)
call ectrans_setup0(kout=nout, kerr=nerr, kverbosity=verbosity, kprgpns=nprgpns, kprgpew=nprgpew, &
  & kprtrw=nprtrw, lduse_mpi=.not.luse_mpi) 
call gstats(1, 1)
!
call gstats(2, 0)
call ectrans_setup(ksmax=nsmax, kdgl=ndgl, kloen=nloen, lduseflt=luseflt) 
call gstats(2, 1)
!
call ectrans_trans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg)


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
call ectrans_print_runtimepars('default')
end if

!===================================================================================================
! Allocate and Initialize spectral arrays
!===================================================================================================

call ectrans_allocate_spectral(nflevl,nspec2,nfld,nsmax)
!!===================================================================================================
!! Allocate gridpoint arrays
!!===================================================================================================
if(.not.allocated(ivset)) allocate(ivset(nflevg))

! Compute spectral distribution
ilev = 0
do jb = 1, nprtrv
  do jlev=1, numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo
!
!
call ectrans_allocate_grid(nproma, ngpblks, nfld, nflevg, lvordiv, luvders, lscders)

!!===================================================================================================
!! Allocate norm arrays
!!===================================================================================================
!
if (lprint_norms .or. ncheck > 0) then
call ectrans_allocate_normdata(nflevl=nflevl, nflevg=nflevg)
call ectrans_calculate_norms(indx=1, nflevl=nflevl, nflevg=nflevg, ivset=ivset,ivsetsc=ivsetsc)
  if (verbosity >= 1) then
    call ectrans_print_norms_init(nout,nflevg)
  endif
endif
!
!!===================================================================================================
!! Setup timers
!!===================================================================================================
!
ztinit = (timef() - ztinit)/1000.0_jprd

if (verbosity >= 0) then
  write(nout,'(" ")')
  write(nout,'(a,i6,a,f9.2,a)') "transform_test initialisation, on",nproc,&
                                & " tasks, took",ztinit," sec"
  write(nout,'(" ")')
endif

if (iters <= 0) call abor1('transform_test:iters <= 0')
call ectrans_allocate_timers(iters)

write(nout,'(a)') '======= Start of spectral transforms  ======='
write(nout,'(" ")')
!
ztloop = timef()
!
!!===================================================================================================
!! Do spectral transform loop
!!===================================================================================================
do jstep = 1, iters
  call gstats(3,0)
! ztstep(jstep) = timef()
  call ectrans_set_ztstep_start(indx=0,jstep=jstep)
!
!  !=================================================================================================
!  ! Do inverse transform
!  !=================================================================================================
!
! ztstep1(jstep) = timef()
  call ectrans_set_ztstep_start(indx=1,jstep=jstep)
  call gstats(4,0)
  if (lvordiv) then
  call ectrans_inv_trans(nproma=nproma,lscders=lscders,luvders=luvders,ivset=ivset,ivsetsc=ivsetsc)
  else
  call ectrans_inv_trans(nproma=nproma,lscders=lscders,ivset=ivset,ivsetsc=ivsetsc)
  endif
  call gstats(4,1)
!
  call ectrans_set_ztstep_end(indx=1,jstep=jstep)
!
!  !=================================================================================================
!  ! While in grid point space, dump the values to disk, for debugging only
!  !=================================================================================================
!
   if (ldump_values) then
!    ! dump a field to a binary file
     call ectrans_dump(jstep, myproc, nproma, ngpblks, nflevg, noutdump)
   endif
 
  !=================================================================================================
  ! Do direct transform
  !=================================================================================================

  call ectrans_set_ztstep_start(indx=2,jstep=jstep)

  call gstats(5,0)
  if (lvordiv) then
  call ectrans_direct_trans(nproma,nfld,ivset,ivsetsc,lvordiv)
  else
  call ectrans_direct_trans(nproma,nfld,ivset,ivsetsc)
  endif
  call gstats(5,1)
  call ectrans_set_ztstep_end(indx=2,jstep=jstep)

  !=================================================================================================
  ! Calculate timings
  !=================================================================================================
  call ectrans_calculate_timings(jstep)
 
   !=================================================================================================
   ! Print norms
   !=================================================================================================
 
   if (lprint_norms) then
     call gstats(6,0)
     call ectrans_calculate_norms(indx=0, nflevl=nflevl, nflevg=nflevg, ivset=ivset,ivsetsc=ivsetsc)
     call ectrans_print_norms_calc(nout, jstep, myproc,nflevg)
    call gstats(6,1)
  else
    call ectrans_print_timestep(nout,jstep)
  endif
  call gstats(3,1)
enddo
!
!!===================================================================================================
!
ztloop = (timef() - ztloop)/1000.0_jprd

write(nout,'(" ")')
write(nout,'(a)') '======= End of spectral transforms  ======='
write(nout,'(" ")')
!
!
if (lprint_norms .or. ncheck > 0) then
     call ectrans_calculate_norms(indx=2, nflevl=nflevl, nflevg=nflevg, ivset=ivset,ivsetsc=ivsetsc)
   if (myproc == 1) then
    call ectrans_print_norms_fin(nout,nflevg,myproc)
   endif
   if (ncheck > 0) then
     ierr = 0
     if (myproc == 1) then
!      ! If the maximum spectral norm error across all fields is greater than 100 times the machine
!      ! epsilon, fail the test
       ierr=ectrans_print_norms_fails(nout,ncheck)
     endif
!
    ! Root rank broadcasts the correctness checker result to the other ranks
    if (luse_mpi) then
      call mpl_broadcast(ierr,kroot=1,ktag=1)
    endif
!
    ! Halt if correctness checker failed
    if (ierr == 1) then
      error stop
     endif
   endif
 endif
!
if (luse_mpi) then
  call ectrans_norms_reduce(ztloop)
endif
call ectrans_compute_time_stats(nproc,iters)
call ectrans_print_time_stats(nout,ztloop,nproc)
!
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
call ectrans_deallocate_grid
!
!!===================================================================================================
!
if (lstats) then
  call gstats(0,1)
  call gstats_print(nout, zaveave, jpmaxstat)
endif

if (lmeminfo) then
  write(nout,*)
  call ec_meminfo(nout, "", mpl_comm, kbarr=1, kiotask=-1, &
      & kcall=1)
endif
enddo !ibackend
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

subroutine get_command_line_arguments(nsmax, cgrid, iters, nfld, nlev, lvordiv, lscders, luvders, &
  &                                   luseflt, nproma, verbosity, ldump_values, lprint_norms, &
  &                                   lmeminfo, nprtrv, nprtrw, ncheck)

  integer, intent(inout) :: nsmax           ! Spectral truncation
  character(len=16), intent(inout) :: cgrid ! Spectral truncation
  integer, intent(inout) :: iters           ! Number of iterations for transform test
  integer, intent(inout) :: nfld            ! Number of scalar fields
  integer, intent(inout) :: nlev            ! Number of vertical levels
  logical, intent(inout) :: lvordiv         ! Also transform vorticity/divergence
  logical, intent(inout) :: lscders         ! Compute scalar derivatives
  logical, intent(inout) :: luvders         ! Compute uv East-West derivatives
  logical, intent(inout) :: luseflt         ! Use fast Legendre transforms
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
      case('--uvders'); luvders = .True.
      case('--flt'); luseflt = .True.
      case('--nproma'); nproma = get_int_value('--nproma', iarg)
      case('--dump-values'); ldump_values = .true.
      case('--norms'); lprint_norms = .true.
      case('--meminfo'); lmeminfo = .true.
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

function cubic_octahedral_gaussian_grid(nsmax) result(cgrid)

  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'O',nsmax+1

end function

!===================================================================================================

subroutine str2int(str, int, stat)

  character(len=*), intent(in) :: str
  integer, intent(out) :: int
  integer, intent(out) :: stat
  read(str, *, iostat=stat) int

end subroutine str2int

subroutine compute_grid_extents()
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
  integer(kind=jpim) :: ib
  integer(kind=jpim) :: jb
  integer(kind=jpim) :: ilev, jlev
  allocate(ivset(nflevg))
  
  ! Compute spectral distribution
  ilev = 0
  do jb = 1, nprtrv
    do jlev=1, numll(jb)
      ilev = ilev + 1
      ivset(ilev) = jb
    enddo
  enddo
  
end subroutine compute_grid_extents
subroutine print_help(unit)

  integer, optional :: unit
  integer :: nout = 6
  if (present(unit)) then
    nout = unit
  endif

  write(nout, "(a)") ""

! if (jprb == jprd) then
!   write(nout, "(a)") "NAME    ectrans-benchmark-dp"
! else
!   write(nout, "(a)") "NAME    ectrans-benchmark-sp"
! end if
! write(nout, "(a)") ""

! write(nout, "(a)") "DESCRIPTION"
! write(nout, "(a)") "        This program tests ecTrans by transforming fields back and forth&
!   & between spectral "
! if (jprb == jprd) then
!   write(nout, "(a)") "        space and grid-point space (double-precision version)"
! else
!   write(nout, "(a)") "        space and grid-point space (single-precision version)"
! end if
! write(nout, "(a)") ""

! write(nout, "(a)") "USAGE"
! if (jprb == jprd) then
!   write(nout, "(a)") "        ectrans-benchmark-dp [options]"
! else
!   write(nout, "(a)") "        ectrans-benchmark-sp [options]"
! end if
! write(nout, "(a)") ""

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
  write(nout, "(a)") ""
  write(nout, "(a)") "DEBUGGING"
  write(nout, "(a)") "    --dump-values       Output gridpoint fields in unformatted binary file"
  write(nout, "(a)") ""

end subroutine print_help

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
subroutine gstats_labels(koffset)
integer koffset        

  call gstats_label(0            , '   ', 'PROGRAM        - Total')
  call gstats_label(1   + koffset, '   ', 'SETUP_TRANS0   - Setup ecTrans')
  call gstats_label(2   + koffset, '   ', 'SETUP_TRANS    - Setup ecTrans handle')
  call gstats_label(3   + koffset, '   ', 'TIME STEP      - Time step')
  call gstats_label(4   + koffset, '   ', 'INV_TRANS      - Inverse transform')
  call gstats_label(5   + koffset, '   ', 'DIR_TRANS      - Direct transform')
  call gstats_label(6   + koffset, '   ', 'NORMS          - Norm comp. (optional)')
  call gstats_label(102 + koffset, '   ', 'LTINV_CTL      - Inv. Legendre transform')
  call gstats_label(103 + koffset, '   ', 'LTDIR_CTL      - Dir. Legendre transform')
  call gstats_label(106 + koffset, '   ', 'FTDIR_CTL      - Dir. Fourier transform')
  call gstats_label(107 + koffset, '   ', 'FTINV_CTL      - Inv. Fourier transform')
  call gstats_label(140 + koffset, '   ', 'SULEG          - Comp. of Leg. poly.')
  call gstats_label(152 + koffset, '   ', 'LTINV_CTL      - M to L transposition')
  call gstats_label(153 + koffset, '   ', 'LTDIR_CTL      - L to M transposition')
  call gstats_label(157 + koffset, '   ', 'FTINV_CTL      - L to G transposition')
  call gstats_label(158 + koffset, '   ', 'FTDIR_CTL      - G to L transposition')
  call gstats_label(400 + koffset, '   ', 'GSTATS         - GSTATS itself')

end subroutine gstats_labels

!===================================================================================================
subroutine ectrans_print_runtimepars(variant)
character(len=*) :: variant
  write(nout,'(" ")')
  write(nout,'(a)')'======= Start of runtime parameters ======='
  write(nout,'("ectrans variant",a)') trim(variant)
  write(nout,'(" ")')
  write(nout,'("nsmax     ",i0)') nsmax
  write(nout,'("grid      ",a)') trim(cgrid)
  write(nout,'("ndgl      ",i0)') ndgl
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
  write(nout,'("nproma    ",i0)') nproma
  write(nout,'("ngpblks   ",i0)') ngpblks
  write(nout,'("nspec2    ",i0)') nspec2
  write(nout,'("nspec2g   ",i0)') nspec2g
  write(nout,'("luseflt   ",l)') luseflt
  write(nout,'("lvordiv   ",l)') lvordiv
  write(nout,'("lscders   ",l)') lscders
  write(nout,'("luvders   ",l)') luvders
  write(nout,'(" ")')
  write(nout,'(a)') '======= End of runtime parameters ======='
  write(nout,'(" ")')
end subroutine ectrans_print_runtimepars
end program transform_test

!===================================================================================================
