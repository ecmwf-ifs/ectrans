! (C) Crown Copyright 2025- Met Office.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

program test_split_mpi_comm

USE PARKIND1, ONLY: JPRM, JPIM
USE MPL_MODULE  ,ONLY : MPL_INIT, MPL_END, MPL_BARRIER, MPL_MYRANK, MPL_NPROC, &
                        MPL_COMM_SPLIT
USE ABORT_TRANS_MOD, ONLY : ABORT_TRANS
USE mpl_data_module, ONLY : MPLUSERCOMM, LMPLUSERCOMM
USE mpl_mpif

implicit none

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"

integer(kind=JPIM), parameter, dimension(2) :: truncations = [79, 188]

! MODE
integer(kind=JPIM), parameter, dimension(2) :: Ms = [1, 3]
integer(kind=JPIM), parameter, dimension(2) :: Ns = [2, 4]


integer(kind=JPIM) :: num_spectral_elements, num_grid_points
integer(kind=JPIM) :: g_num_spectral_elements, g_num_grid_points  ! global
integer(kind=JPIM) :: local_spectral_coefficient_index
integer(kind=JPIM) :: ierror
integer(kind=JPIM) :: i
integer(kind=JPIM) :: split_num_ranks, split_rank
integer(kind=JPIM) :: world_num_ranks, world_rank
integer(kind=JPIM) :: split_colour, split_key
integer(kind=JPIM) :: split_comm
integer(kind=JPIM) :: truncation
integer(kind=JPIM) :: M, N
integer(kind=JPIM) :: num_latitudes, num_longitudes

! Number of grid points on each rank
integer(kind=JPIM) :: grid_partition_size_local(1)
integer(kind=JPIM), allocatable :: displs(:)
integer(kind=JPIM), allocatable :: grid_partition_sizes(:)

integer(kind=JPIM), allocatable :: spectral_indices(:)

real(kind=JPRM), allocatable :: spectral_field(:,:)
real(kind=JPRM), allocatable :: grid_point_field(:,:,:)

! NOTE: 1 Dimensional as this is a field simply used to write to file output.
real(kind=JPRM), allocatable :: g_grid_point_field(:)

character(len=1024) :: filename


call MPI_Init(ierror)
call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierror)

split_colour = get_split_group()
split_key = world_rank
call MPI_Comm_split(MPI_COMM_WORLD, split_colour, split_key, split_comm, ierror)

print*,"=== Rank ", world_rank, ", Setup on group", split_colour, "==="

! Set MPL comm. Global variables from fiat MPL_DATA_MODULE.
! Acts like arguments to MPL_INIT.
LMPLUSERCOMM = .true.
MPLUSERCOMM = split_comm
call MPL_INIT()

split_rank = MPL_MYRANK()
split_num_ranks = MPL_NPROC()

! Assert that the split comm is smaller than WORLD.
if (split_num_ranks >= world_num_ranks) then
  call ABORT_TRANS("ERROR: Split communicator not smaller than MPI_COMM_WORLD.")
end if

print*,"=== Local rank ", split_rank, ", on group", split_colour, "size", split_num_ranks, "==="

call setup_trans0(KPRINTLEV=0, LDMPOFF=.false., &
!                 Split grid NS            Split spectral
                  KPRGPNS=split_num_ranks, KPRTRW=split_num_ranks)

! Different transform based on the colour.
truncation = truncations(split_colour + 1)

num_latitudes = 2*(truncation + 1)
num_longitudes = num_latitudes*2

if (split_rank == 1) print*,"Colour ", split_colour, &
  " GLOBAL NUM LON =>", num_longitudes, "| LAT =>", num_latitudes

call setup_trans(KSMAX=truncation, KDGL=num_latitudes)

! Get function space sizes
call trans_inq(KSPEC2=num_spectral_elements, KGPTOT=num_grid_points)
call trans_inq(KSPEC2G=g_num_spectral_elements, KGPTOTG=g_num_grid_points)

allocate(spectral_field(1, num_spectral_elements))
allocate(grid_point_field(num_grid_points, 1, 1))

! Get spectral indices
allocate(spectral_indices(0:truncation))
call trans_inq(KASM0=spectral_indices)

! select mode
M = Ms(split_colour + 1)
N = Ns(split_colour + 1)
local_spectral_coefficient_index = spectral_indices(M) + 2*(N - M) + 1

spectral_field(:,:) = 0.0

if (local_spectral_coefficient_index > 0) then
  spectral_field(1,local_spectral_coefficient_index) = 1.0
end if

call inv_trans(PSPSCALAR=spectral_field, PGP=grid_point_field)

! -------------- Gather the result on the root (0) and write to file -----------------

! Get counts from each PE.
grid_partition_size_local(1) = num_grid_points
allocate(grid_partition_sizes(split_num_ranks))
grid_partition_sizes = 0

call MPI_Gather(grid_partition_size_local, 1, MPI_INT, &
                grid_partition_sizes, 1, MPI_INT, &
                0, split_comm, ierror)
if (ierror /= 0) then
  print*,"MPI ERROR => ", ierror
  call ABORT_TRANS("MPI ERROR")
end if

if (split_rank == 1) then
  ! Allocate a global field.
  allocate(g_grid_point_field(g_num_grid_points))

  ! Make displacement arrays
  allocate(displs(split_num_ranks))
  displs = 0
  do i=2, split_num_ranks
    displs(i) = displs(i - 1) + grid_partition_sizes(i - 1)
  end do
end if

call MPI_Gatherv(grid_point_field(:,1,1), num_grid_points, MPI_FLOAT, &
                 g_grid_point_field, grid_partition_sizes, displs, MPI_FLOAT, &
                 0, split_comm, ierror)
if (ierror /= 0) then
  print*,"Gatherv MPI ERROR => ", ierror
  call ABORT_TRANS("MPI ERROR")
end if

if (split_rank == 1) then
  write(filename, "(A22,I0,A4)") "grid_point_field_trunc_", truncation, ".dat"
  open(7, file=filename, form="unformatted")
  write(7) g_grid_point_field(:)
  close(7)

  print*,"Colour", split_colour, "finished and written to file: "//trim(filename)
end if

call MPL_END()
call MPI_Finalize()

CONTAINS

! Get the colour of comm for this rank.
function get_split_group() result(group)
  implicit none
  ! return
  integer(kind=JPIM) :: group

  integer(kind=JPIM) :: rank, world_size, ierror
  real(kind=JPRM) :: rank_ratio

  call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

  ! ----------------------------------------------
  ! Uneven splitting.
  ! ----------------------------------------------
  rank_ratio = real(rank + 1, kind=JPRM) / real(world_size, kind=JPRM)

  ! Split X%
  if (rank_ratio <= 0.25_jprm) then
    group = 0
  else
    group = 1
  end if

end function get_split_group

end program
