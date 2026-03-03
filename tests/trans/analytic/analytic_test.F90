! (C) Copyright 2023- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program analytic_test

use parkind1, only: jpim, jprb, jprd, jprm
use mpl_module
! use analytic_solutions_mod, only: analytic_init, analytic_end, &
! & buffer_legendre_polynomials_supolf, &
! & buffer_legendre_polynomials_ectrans, check_legendre_polynomials, &
! & compute_analytic_solution, compute_analytic_eastwest_derivative, &
! & compute_analytic_northsouth_derivative, gelam, gelat, init_check_fields, &
! & close_check_fields, check_gp_fields, check_sp_fields, compute_analytic_uv, &
! & compute_analytic_uv_derivative_ew

implicit none

! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT

! Default parameters
integer(kind=jpim) :: nsmax   = 21  ! Spectral truncation
integer(kind=jpim) :: nfld    = 1   ! Number of scalar fields
logical            :: limag = .false. ! test imaginary part of spectral data
integer(kind=jpim) :: ndgl ! Number of latitudes
integer(kind=jpim) :: nspec2
integer(kind=jpim) :: ngptot
integer(kind=jpim) :: ngptotg
integer(kind=jpim) :: nspec2g
integer(kind=jpim) :: ja
integer(kind=jpim) :: ib
integer(kind=jpim) :: jprtrv
integer(kind=jpim) :: n_regions_ns
integer(kind=jpim) :: n_regions_ew

integer(kind=jpim), allocatable :: nloen(:)
integer(kind=jpim) :: myproc

! Spectral and grid point arrays
real(kind=jprb), allocatable :: zspscalar(:,:)
real(kind=jprb), allocatable :: zgp(:,:,:)

! Spectral space data structures
real(kind=jprd), allocatable :: zsph_analytic(:,:)

logical :: luseflt = .false. ! Use fast legendre transforms
real(kind=jprd) :: rtolerance = 1e-9 ! maximum relative lmax error tolerance for
                                     ! passing analytyic solution tests in double precision
                                     ! Value for single precision is set at the beginning
                                     ! of the execution.

! Verbosity level (-1, 0 or 1)
integer :: verbosity = -1

integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nprgpns ! Grid-point decomp
integer(kind=jpim) :: nprgpew ! Grid-point decomp
integer(kind=jpim) :: nprtrv = 0 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
integer(kind=jpim) :: mysetv

integer(kind=jpim), allocatable :: num_local_levels_all(:), ivset(:)

integer(kind=jpim) :: nfld_local, num_rest
integer(kind=jpim) :: i
integer(kind=jpim) :: isqr
integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: ngpblks
integer(kind=jpim) :: ilev, jlev
integer(kind=jpim) :: jsetv
integer(kind=jpim) :: m, n, imag_idx
integer, external :: ec_mpirank
logical :: luse_mpi = .true.
character(len=16) :: cgrid = ''

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "abor1.intfb.h"

!===================================================================================================

luse_mpi = detect_mpirun()
if (jprb == jprm) rtolerance = 1e-3 ! tolerance for single precision
! Setup
call get_command_line_arguments(nsmax, cgrid, nfld, luseflt, nproma, verbosity, nprtrv, nprtrw, limag, rtolerance)
if (cgrid == '') cgrid = cubic_full_grid(nsmax)
call parse_grid(cgrid, ndgl, nloen)

write(nout,'(a)')'======= ecTrans analytic test ======='
write(nout,'(a,i0)') 'Spectral truncation: ', nsmax
write(nout,'(a,a)') 'Grid: ', trim(cgrid)
write(nout,'(a,i0)') 'Number of scalar fields: ', nfld

!===================================================================================================

if (luse_mpi) then
  call mpl_init(ldinfo=(verbosity >= 1))
  nproc  = mpl_nproc()
  myproc = mpl_myrank()
else
  nproc = 1
  myproc = 1
endif

!===================================================================================================

! only output to stdout on pe 1
if (myproc /= 1) then
  open(unit=nout, file='/dev/null')
endif

!===================================================================================================

! Compute nprgpns and nprgpew
! This version selects most square-like distribution
isqr = int(sqrt(real(nproc, jprb)))
do ja = isqr, nproc
  ib = nproc / ja
  if (ja * ib == nproc) then
    nprgpns = max(ja, ib)
    nprgpew = min(ja, ib)
    exit
  endif
enddo

! Compute nprtrv and nprtrw if not provided on the command line
if (nprtrv > 0 .or. nprtrw > 0) then
  if (nprtrv == 0) nprtrv = nproc / nprtrw
  if (nprtrw == 0) nprtrw = nproc / nprtrv
  if (nprtrw * nprtrv /= nproc) call abor1('analytic_test: nprtrw * nprtrv /= nproc')
else
  do jprtrv = 4, nproc
    nprtrv = jprtrv
    nprtrw = nproc / nprtrv
    if (nprtrv * nprtrw /= nproc) cycle
    if (nprtrv > nprtrw) exit
  enddo
  ! Go for approx square partition for backup
  if (nprtrv * nprtrw /= nproc .or. nprtrv > nprtrw) then
    isqr = int(sqrt(real(nproc, jprb)))
    do ja = isqr, nproc
      ib = nproc / ja
      if (ja * ib == nproc) then
        nprtrw = max(ja, ib)
        nprtrv = min(ja, ib)
        exit
      endif
    enddo
  endif
endif

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

if (verbosity >= 1) write(nout,'(a)')'======= Setup ecTrans ======='

call setup_trans0(kout=nout, kerr=nerr, kprintlev=merge(2, 0, verbosity == 1), kprgpns=nprgpns, kprgpew=nprgpew, &
  &               kprtrw=nprtrw, ldalloperm=.true.,                      &
  &               ldmpoff=.not.luse_mpi, k_regions_ns=n_regions_ns, k_regions_ew=n_regions_ew)

call setup_trans(ksmax=nsmax, kdgl=ndgl, kloen=nloen, lduseflt=luseflt)

call trans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg)

if (nproma == 0) then ! no blocking (default when not specified)
  nproma = ngptot
endif

! Calculate number of NPROMA blocks
ngpblks = (ngptot - 1)/nproma+1

if (nproc == 1) then
  mysetv = 1
else
  call trans_inq(kmysetv=mysetv)
endif

! Determine number of local levels for fourier and legendre calculations
! based on the values of nfld and nprtrv
allocate(num_local_levels_all(nprtrv))
nfld_local = nfld / nprtrv ! Integer division
num_rest = nfld - nfld_local * nprtrv ! Tail block
do jsetv = 1, nprtrv
  num_local_levels_all(jsetv) = nfld_local
  if (jsetv == nprtrv) then
    num_local_levels_all(jsetv) = num_local_levels_all(jsetv) + num_rest
  endif
enddo

nfld_local = num_local_levels_all(mysetv)

!===================================================================================================
! Allocate arrays
!===================================================================================================

allocate(ivset(nfld))

! Compute spectral distribution
ilev = 0
do jsetv = 1, nprtrv
  do jlev = 1, num_local_levels_all(jsetv)
    ilev = ilev + 1
    ivset(ilev) = jsetv
  enddo
enddo

allocate(zspscalar(nfld_local,nspec2))
allocate(zgp(nproma,nfld,ngpblks))

! Allocate arrays for analytic solutions
allocate(zsph_analytic(nproma,ngpblks))

! Compute geographic longitude gelam and latitude gelat:
! call analytic_init(nproma, ngpblks, ndgl, n_regions_ns, n_regions_ew, nloen)
! call buffer_legendre_polynomials_supolf(nsmax)

!===================================================================================================
! Perform tests
!===================================================================================================

! Loop over all wavenumbers (check actually tested wavenumber inside)
do n = 0, nsmax
  do m = 0, n
    do imag_idx = 0, 1
      call initialize_spectral_array(nsmax, m, n, zspscalar)

      !=================================================================================================
      ! Compute analytic solutions
      !=================================================================================================

      ! call compute_analytic_solution(nproma, ngpblks, nsmax, ngptot, m, n, li, zsph_analytic)

      !=================================================================================================
      ! Do inverse transform
      !=================================================================================================

      call inv_trans(pspscalar=zspscalar, kproma=nproma, kvsetsc=ivset, pgp=zgp)

      !=================================================================================================
      ! Do direct transform
      !=================================================================================================

      call dir_trans(pgp=zgp, kproma=nproma, kvsetsc=ivset, pspscalar=zspscalar)
    end do
  end do
end do

!===================================================================================================
! Cleanup
!===================================================================================================

deallocate(zsph_analytic, zgp, zspscalar)

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

function get_real_value(cname, iarg) result(value)

  real :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg
  character(len=128) :: carg
  integer :: stat

  carg = get_str_value(cname, iarg)
  call str2real(carg, value, stat)

  if (stat /= 0) then
    call parsing_failed("Invalid argument for " // trim(cname) // ": " // trim(carg))
  end if

end function

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

subroutine get_command_line_arguments(nsmax, cgrid, nfld, luseflt, nproma, verbosity, nprtrv, nprtrw, limag, rtolerance)

  integer, intent(inout) :: nsmax           ! Spectral truncation
  character(len=16), intent(inout) :: cgrid ! Grid
  integer, intent(inout) :: nfld            ! Number of scalar fields
  logical, intent(inout) :: luseflt         ! Use fast Legendre transforms
  integer, intent(inout) :: nproma          ! NPROMA
  integer, intent(inout) :: verbosity       ! Level of verbosity
  integer, intent(inout) :: nprtrv          ! Size of V set (spectral decomposition)
  integer, intent(inout) :: nprtrw          ! Size of W set (spectral decomposition)
  logical, intent(inout) :: limag           ! test imaginary part
  real(jprd), intent(inout) :: rtolerance      ! relative error tolerance for analytic solutions

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
      ! Parse spectral truncation argument
      case('-t', '--truncation')
        nsmax = get_int_value('-t', iarg)
        if (nsmax < 1) then
          call parsing_failed("Invalid argument for -t: must be > 0")
        end if
      case('-g', '--grid'); cgrid = get_str_value('-g', iarg)
      case('-f', '--nfld'); nfld = get_int_value('-f', iarg)
      case('--imaginary'); limag = .true.
      case('--flt'); luseflt = .True.
      case('--nproma'); nproma = get_int_value('--nproma', iarg)
      case('--nprtrv'); nprtrv = get_int_value('--nprtrv', iarg)
      case('--nprtrw'); nprtrw = get_int_value('--nprtrw', iarg)
      case('--tolerance'); rtolerance = get_real_value('--tolerance', iarg)
      case default
        call parsing_failed("Unrecognised argument: " // trim(carg))

    end select
    iarg = iarg + 1
  end do

end subroutine get_command_line_arguments

!===================================================================================================

function cubic_octahedral_gaussian_grid(nsmax) result(cgrid)

  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'O',nsmax+1

end function

!===================================================================================================

function cubic_full_grid(nsmax) result(cgrid)

  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'F',nsmax+1

end function

!===================================================================================================

subroutine str2int(str, int, stat)

  character(len=*), intent(in) :: str
  integer, intent(out) :: int
  integer, intent(out) :: stat
  read(str, *, iostat=stat) int

end subroutine str2int

!===================================================================================================

subroutine str2real(str, real, stat)

  character(len=*), intent(in) :: str
  real, intent(out) :: real
  integer, intent(out) :: stat
  read(str, *, iostat=stat) real

end subroutine str2real

!===================================================================================================

subroutine print_help(unit)

  integer, optional :: unit
  integer :: nout = 6
  if (present(unit)) then
    nout = unit
  endif

  write(nout, "(a)") ""

  if (jprb == jprd) then
    write(nout, "(a)") "NAME    ectrans-benchmark-dp"
  else
    write(nout, "(a)") "NAME    ectrans-benchmark-sp"
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
    write(nout, "(a)") "        ectrans-benchmark-dp [options]"
  else
    write(nout, "(a)") "        ectrans-benchmark-sp [options]"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "OPTIONS"
  write(nout, "(a)") "    -h, --help          Print this message"
  write(nout, "(a)") "    -v                  Run with verbose output"
  write(nout, "(a)") "    -t, --truncation T  Run with this triangular spectral truncation"
  write(nout, "(a)") "                        (default = 79)"
  write(nout, "(a)") "    -g, --grid GRID     Run with this grid. Possible values: O<N>, F<N>"
  write(nout, "(a)") "                        If not specified, O<N> is used with N=truncation+1"
  write(nout, "(a)") "                        (cubic relation)"
  write(nout, "(a)") "    -n, --niter NITER   Run for this many inverse/direct transform"
  write(nout, "(a)") "                        iterations (default = 10)"
  write(nout, "(a)") "    -f, --nfld NFLD     Number of scalar fields (default = 1)"
  write(nout, "(a)") "    --imaginary         Test imaginary part"
  write(nout, "(a)") "    --tolerance         Test is passed if largest relative lmax-error is"
  write(nout, "(a)") "                        smaller than this tolerance (real value)"
  write(nout, "(a)") "    --vordiv            Also transform vorticity-divergence to wind"
  write(nout, "(a)") "    --scders            Compute scalar derivatives (default off)"
  write(nout, "(a)") "    --flt               Run with fast Legendre transforms (default off)"
  write(nout, "(a)") "    --nproma NPROMA     Run with NPROMA (default no blocking: NPROMA=ngptot)"
  write(nout, "(a)") "    --norms             Calculate and print spectral norms of transformed"
  write(nout, "(a)") "                        fields"
  write(nout, "(a)") "                        The computation of spectral norms will skew overall"
  write(nout, "(a)") "                        timings"
  write(nout, "(a)") "    --nprtrv            Size of V set in spectral decomposition"
  write(nout, "(a)") "    --nprtrw            Size of W set in spectral decomposition"
  write(nout, "(a)") "    -c, --check VALUE   The multiplier of the machine epsilon used as a"
  write(nout, "(a)") "                        tolerance for correctness checking"
  write(nout, "(a)") ""

end subroutine print_help

!===================================================================================================

subroutine initialize_spectral_array(nsmax, kzonal, ktotal, pspscalar)

  integer,            intent(in)  :: nsmax          ! Spectral truncation
  integer,            intent(in)  :: kzonal         ! Zonal wavenumber
  integer,            intent(in)  :: ktotal         ! Total wavenumber
  real(kind=jprb),    intent(out) :: pspscalar(:,:) ! Input spectral array

  integer :: jfld
  integer :: index, num_my_zon_wns
  integer, allocatable :: my_zon_wns(:), nasm0(:)

  ! Get zonal wavenumbers this rank is responsible for
  call trans_inq(knump=num_my_zon_wns)
  allocate(my_zon_wns(num_my_zon_wns))
  call trans_inq(kmyms=my_zon_wns)

  ! First initialise all spectral coefficients to zero
  pspscalar(:,:) = 0.0

  ! If rank is responsible for the chosen zonal wavenumber...
  if (any(my_zon_wns == kzonal)) then
    ! Get array of spectral array addresses (this maps (m, n=m) to array index)
    allocate(nasm0(0:nsmax))
    call trans_inq(kasm0=nasm0)

    ! Find out local array index of chosen spherical harmonic
    index = nasm0(kzonal) + 2 * (ktotal - kzonal)

    ! Set just that element to a constant value
    do jfld = 1, size(pspscalar, 1)
      pspscalar(jfld, index) = 1.0_jprb
    end do
  end if

end subroutine initialize_spectral_array

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

end program analytic_test

!===================================================================================================
