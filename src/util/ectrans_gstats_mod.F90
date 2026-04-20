module ectrans_gstats_mod

implicit none
!private
public :: ectrans_gstats_new_region
public :: ectrans_gstats_labels
public :: ectrans_gstats_init
public :: ectrans_gstats_end
public :: ectrans_gstats_enable
public :: ectrans_gstats_configuration

logical :: lstats_ = .false.
integer, save :: gstats_counter = 10

type :: ectrans_gstats_configuration
  logical :: lstats = .false.
  logical :: lstatscpu = .false.
  logical :: lsyncstats = .false.
  logical :: ldetailed_stats = .false.
  logical :: lbarrier_stats = .false.
  logical :: lbarrier_stats2 = .false.
  logical :: lstats_omp = .false.
  logical :: lstats_comms = .false.
  logical :: lstats_mem = .false.
  integer :: nstats_mem = 0
  logical :: lstats_alloc = .false.
  logical :: ltrace_stats = .false.
  integer :: ntrace_stats = 0
  integer :: nprnt_stats = 1
  logical :: lxml_stats = .false.
end type

contains

function new_gstats_id() result(gstats_id)
  integer :: gstats_id
  gstats_id = gstats_counter
  gstats_counter = gstats_counter + 1
end function

function ectrans_gstats_new_region(category,description) result(gstats_id)
  character(len=*), optional, intent(in) :: category
  character(len=*), optional, intent(in) :: description
  character(len=3)   :: category_local
  character(len=256) :: description_local
  integer :: gstats_id
  if (.not. lstats_) then
    gstats_id = -1
    return
  endif
  gstats_id = new_gstats_id()
  if (present(category)) then
    category_local = category(1:max(3,len_trim(category)))
  else
    category_local = ''
  end if
  if (present(description)) then
    description_local = description(1:min(len(description_local),len_trim(description)))
  else
    description_local = ''
  end if
  write(0,*) 'Registering GSTATS region ', gstats_id, ': [', category_local, '] ', trim(description_local)
  call gstats_label(gstats_id, trim(category_local), trim(description_local))
end function

! Assign GSTATS labels to the main regions of ecTrans
subroutine ectrans_gstats_labels
  
  ! From gstats_label_ifs.F90

  !   counters 101 to 200 :  trans package
  call gstats_label(102,'TRS','LTINV_CTL   - INVERSE LEGENDRE TRANSFORM')
  call gstats_label(103,'TRS','LTDIR_CTL   - DIRECT LEGENDRE TRANSFORM')
  call gstats_label(104,'TRS','LTINV_CTLAD - ADJ. INVERSE LEGENDRE TRANSFORM')
  call gstats_label(105,'TRS','LTDIR_CTLAD - ADJ. DIRECT LEGENDRE TRANSFORM')
  call gstats_label(106,'TRS','FTDIR_CTL   - DIRECT FOURIER TRANSFORM')
  call gstats_label(107,'TRS','FTINV_CTL   - INVERSE FOURIER TRANSFORM')
  call gstats_label(108,'TRS','OMP in GPC  - TEST_ADJOINT ')
  call gstats_label(110,'TRS','OMP in SCALPRODSP  - TEST_ADJOINT ')
  call gstats_label(132,'TRS','FTINV_CTLAD - ADJ. INVERSE FOURIER TRANSFORM')
  call gstats_label(133,'TRS','FTDIR_CTLAD - ADJ. DIRECT FOURIER TRANSFORM')
  call gstats_label(135,'TRS','OMP in GPCAD - TEST_ADJOINT ')
  call gstats_label(140,'TRS','SULEG       - COMP. OF LEGENDRE POL.')
  call gstats_label(152,'TRS','LTINV_CTL   - M TO L TRANSPOSITION')
  call gstats_label(153,'TRS','LTDIR_CTL   - L TO M TRANSPOSITION')
  call gstats_label(157,'TRS','FTINV_CTL   - L TO G TRANSPOSITION')
  call gstats_label(158,'TRS','FTDIR_CTL   - G TO L TRANSPOSITION')
  call gstats_label(180,'TRS','LTINV_CTLAD - L TO M TRANSPOSITION')
  call gstats_label(181,'TRS','LTDIR_CTLAD - M TO L TRANSPOSITION')
  call gstats_label(182,'TRS','FTINV_CTLAD - G TO L TRANSPOSITION')
  call gstats_label(183,'TRS','FTDIR_CTLAD - L TO G TRANSPOSITION')
  call gstats_label(190,'TRS','SUTRLE      - COMMUNICATE LEG.POL.')

  !   counters 400-401 are reserved for gstats itself
  call gstats_label(400,'   ','GSTATS         - GSTATS itself')
  call gstats_label(401,'   ','GSTATS HOOK')

  !   It looks like also counters were added in GPU code path in range 410 to 450
  !   TODO: assign labels to these counters as well, or consider renumbering to have
  !   a more systematic assignment of counter ranges to different parts of the code

  ! barrier counter 700 ---> 800

  call gstats_label(761,'GBR','GBAR IN TRGTOL          ')
  call gstats_label(762,'GBR','GBAR IN TRLTOG          ')
  call gstats_label(763,'GBR','GBAR IN TRLTOM          ')
  call gstats_label(764,'GBR','GBAR IN TRMTOL          ')
  call gstats_label(783,'BAR','BARRIER IN SUTRLE')
  call gstats_label(784,'BAR','BARRIER IN GATH_GRID_CTL')
  call gstats_label(785,'BAR','BARRIER IN GATH_SPEC_CONTROL')
  call gstats_label(786,'BAR','BARRIER IN DIST_GRID_CTL')
  call gstats_label(787,'BAR','BARRIER IN DIST_SPEC_CONTROL')
  call gstats_label(788,'GBR','GBAR IN GATH_SPEC_CONTROL')
  call gstats_label(789,'GBR','GBAR IN GATH_GRID_CTL')
  call gstats_label(790,'GBR','GBAR IN DIST_SPEC_CONTROL')
  call gstats_label(791,'GBR','GBAR IN DIST_GRID_CTL')
  call gstats_label(795,'GBR','BARRIER IN SUSTAONL')
  call gstats_label(798,'BAR','BARRIER IN SULEG')

  !   counters 800 to 900 :  trans package
  call gstats_label(801,'MPL','MPI IN SUTRLE_MOD ')
  call gstats_label(803,'MPL','TRGTOL_COMMS')
  call gstats_label(804,'MPL','TRGTOL_COMMS (GPNORM)')
  call gstats_label(805,'MPL','TRLTOG_COMMS')
  call gstats_label(806,'MPL','TRLTOM_COMMS')
  call gstats_label(807,'MPL','TRMTOL_COMMS')
  call gstats_label(809,'MPL','GATH_GRID_CTL_COMMS')
  call gstats_label(810,'MPL','GATH_SPEC_CONTROL_COMMS')
  call gstats_label(811,'MPL','DIST_GRID_CTL_COMMS')
  call gstats_label(812,'MPL','DIST_SPEC_CONTROL_COMMS')
  call gstats_label(813,'MPL','EVCOST')
  call gstats_label(814,'MPL','TRANS SUSTAONL')
  call gstats_label(815,'MPL','GPNORM_TRANS')
  call gstats_label(816,'MPL','GPNORM_TRANS')
  call gstats_label(817,'MPL','SUOBSCOR')
  call gstats_label(851,'MPL','SULEG - SUPOLF')
  call gstats_label(852,'MPL','SULEG - CONSTRUCT_BUTTERFLY')

  call gstats_label(1141,'OMP','SUTRLE       ')
  call gstats_label(1251,'OMP','SULEG - SUPOLF')
  call gstats_label(1252,'OMP','SULEG - CONSTRUCT_BUTTERFLY')
  call gstats_label(1253,'OMP','CONSTRUCT BUTTERFLY ')

  call gstats_label(1429,'OMP','GPNORM_TRANS ')

  !   counters 1600 to 1700 :  trans package  for OMP
  call gstats_label(1601,'OMP','TRGTOL LOCAL ')
  call gstats_label(1602,'OMP','TRGTOL PACK  ')
  call gstats_label(1603,'OMP','TRGTOL UNPACK')
  call gstats_label(1604,'OMP','TRLTOG LOCAL ')
  call gstats_label(1605,'OMP','TRLTOG PACK  ')
  call gstats_label(1606,'OMP','TRLTOG UNPACK')
  call gstats_label(1607,'OMP','TRLTOM       ')
  call gstats_label(1608,'OMP','TRMTOL       ')
  call gstats_label(1639,'OMP','FTINV_CTL    ')
  call gstats_label(1640,'OMP','FTDIR_CTL    ')
  call gstats_label(1641,'OMP','FTINV_CTLAD  ')
  call gstats_label(1642,'OMP','FTDIR_CTLAD  ')
  call gstats_label(1643,'OMP','GATH_GRID_CTL ')
  call gstats_label(1644,'OMP','GATH_SPEC_CONTROL ')
  call gstats_label(1645,'OMP','LTDIR_CTL   - DIRECT LEGENDRE TRANSFORM')
  call gstats_label(1646,'OMP','LTDIR_CTLAD - ADJ. DIRECT LEGENDRE TRANSFORM')
  call gstats_label(1647,'OMP','LTINV_CTL   - INVERSE LEGENDRE TRANSFORM')
  call gstats_label(1648,'OMP','LTINV_CTLAD - ADJ. INVERSE LEGENDRE TRANSFORM')
  call gstats_label(1650,'OMP','SUGAW_MOD    ')
  call gstats_label(1651,'OMP','SPNORMD      ')
  call gstats_label(1655,'OMP','SPCSI LDONEM=F SIGAM')
  call gstats_label(1656,'OMP','SPCSI LDONEM=F')
  call gstats_label(1657,'OMP','SPCSI LDONEM=F SITNU')
  call gstats_label(1658,'OMP','SPCSIAD LDONEM=F SIGAMAD')
  call gstats_label(1659,'OMP','SPCSIAD LDONEM=F')
  call gstats_label(1660,'OMP','SPCSI LDONEM=F MXMAOP')
  call gstats_label(1661,'OMP','WVCOUPLE Triple loop')
  call gstats_label(1662,'OMP','SPCSIAD LDONEM=F SITNUAD')
  call gstats_label(1663,'OMP','DIST_GRID_CTL ')
  call gstats_label(1664,'OMP','SPCSIAD LDONEM=F MXMAOP')

  !   counters 1800 to 1810 :  trans package  for  serial part
  call gstats_label(1801,'SER','SULEG ')
  call gstats_label(1802,'SER','SETUP_TRANS ')
  call gstats_label(1803,'SER','SUTRLE ')
  call gstats_label(1804,'SER','DIST_SPEC_CONTROL_SERIAL')
  call gstats_label(1805,'SER','TRGTOL init')
  call gstats_label(1806,'SER','TRLTOG init')
  call gstats_label(1807,'SER','INV_TRANS init')
  call gstats_label(1808,'SER','DIR_TRANS init')
  call gstats_label(1809,'SER','INV_TRANSAD init')
  call gstats_label(1810,'SER','DIR_TRANSAD init')

  call gstats_label(1901,'SER','CONSTRUCT BUTTERFLY ')

end subroutine ectrans_gstats_labels


subroutine ectrans_gstats_init(config)
  use ectrans_mpi_mod, only : ectrans_mpi_world_rank, ectrans_mpi_world_size
  use ec_parkind, only : jpim
  implicit none
  type(ectrans_gstats_configuration), intent(in) :: config

  integer(kind=jpim) :: nproc, myproc, jproc
  integer(kind=jpim), allocatable :: nprcids(:)
  logical :: lstats_omp
  logical :: lstats_comms
  logical :: lstatscpu
  integer(kind=jpim) :: nprnt_stats

#include "gstats_setup.intfb.h"

  lstats_ = config%lstats

  if (lstats_) then
    nproc  = ectrans_mpi_world_size()
    myproc = ectrans_mpi_world_rank() + 1
    if (nproc == 0) then
      call abor1("ectrans_gstats_init: nproc must be > 0 when lstats is true")
    endif

    allocate(nprcids(nproc))
    nprcids = [(jproc, jproc=1,nproc)]

    lstats_omp    = config%lstats_omp
    lstats_comms  = config%lstats_comms
    lstatscpu     = config%lstatscpu
    nprnt_stats   = config%nprnt_stats
    if (config%ldetailed_stats) then
      lstats_omp    = .true.
      lstats_comms  = .true.
      lstatscpu     = .true.
      nprnt_stats   = nproc
    endif

    call gstats_setup(nproc, myproc, nprcids, lstats_,                                    &
      & lstatscpu, config%lsyncstats, config%ldetailed_stats, config%lbarrier_stats,      &
      & config%lbarrier_stats2, lstats_omp, lstats_comms, config%lstats_mem,              &
      & config%nstats_mem, config%lstats_alloc, config%ltrace_stats, config%ntrace_stats, &
      & config%nprnt_stats, config%lxml_stats)
    call gstats_psut

    call gstats(0,0)
    call gstats_label(0,   '   ', 'PROGRAM        - Total')
  end if

  call ectrans_gstats_labels
end subroutine

subroutine ectrans_gstats_print()
  use ectrans_log_mod, only : nout
  use yomgstats, only: jpmaxstat
  use ec_parkind, only : jprd
  implicit none
  real(kind=jprd) :: zaveave(0:jpmaxstat)
  if (lstats_) then
    call gstats_print(nout, zaveave, jpmaxstat)
  endif
end subroutine

subroutine ectrans_gstats_end()
  implicit none
  if (lstats_) then
    call gstats(0,1)
    call ectrans_gstats_print()
  endif
end subroutine

subroutine ectrans_gstats_enable(lenable)
  use yomgstats, only: lstats
  logical, intent(in), optional :: lenable
  if (lstats_) then
    lstats = .true.
    if (present(lenable)) then
      lstats = lenable
    endif
  endif
end subroutine

end module