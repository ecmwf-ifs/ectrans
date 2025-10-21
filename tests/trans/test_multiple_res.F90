! (C) Crown Copyright 2025- Met Office.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!

program test_multiple_res

USE PARKIND1, ONLY: JPRM, JPIM
USE MPL_MODULE  ,ONLY : MPL_INIT, MPL_END, MPL_BARRIER, MPL_MYRANK, MPL_NPROC, &
                        MPL_COMM_SPLIT
USE ABORT_TRANS_MOD, ONLY : ABORT_TRANS
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
integer(kind=JPIM) :: mode_index
integer(kind=JPIM) :: ierror
integer(kind=JPIM) :: i
integer(kind=JPIM) :: res_idx
integer(kind=JPIM) :: num_ranks, rank
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

real(kind=JPRM), allocatable :: g_grid_point_field(:)

character(len=1024) :: filename

integer(kind=JPIM), dimension(2) :: handles


call MPL_INIT()
num_ranks = MPL_NPROC()
rank = MPL_MYRANK()
print*, "=== Local rank ", rank, "size", num_ranks, "==="

!                 Num resolutions                             Split grid NS      Split spectral
call setup_trans0(KMAX_RESOL=2, KPRINTLEV=0, LDMPOFF=.false., KPRGPNS=num_ranks, KPRTRW=num_ranks)

do res_idx = 1, 2
  truncation = truncations(res_idx)

  num_latitudes = 2*(truncation + 1)
  num_longitudes = num_latitudes*2
  print*, ">>> GLOBAL NUM LON =>", num_longitudes, "| LAT =>", num_latitudes

  !                 OUT
  call setup_trans(KRESOL=handles(res_idx), KSMAX=truncation, KDGL=num_latitudes)

  ! Get function space sizes
  call trans_inq(KRESOL=handles(res_idx), KSPEC2=num_spectral_elements, KGPTOT=num_grid_points, &
  ! Global function space sizes.
                                          KSPEC2G=g_num_spectral_elements, KGPTOTG=g_num_grid_points)
  print*,"Num spec = ", num_spectral_elements, "| Num grid points = ", num_grid_points, g_num_grid_points

  if (allocated(spectral_field)) deallocate(spectral_field)
  allocate(spectral_field(1, num_spectral_elements))
  if (allocated(grid_point_field)) deallocate(grid_point_field)
  allocate(grid_point_field(num_grid_points, 1, 1))

  ! Get spectral indices
  if (allocated(spectral_indices)) deallocate(spectral_indices)
  allocate(spectral_indices(0:truncation))
  call trans_inq(KRESOL=handles(res_idx), KASM0=spectral_indices)

  ! select mode
  M = Ms(res_idx)
  N = Ns(res_idx)
  mode_index = spectral_indices(M) + 2*(N - M) + 1

  spectral_field(:,:) = 0.0

  if (mode_index > 0) then
    spectral_field(1,mode_index) = 1.0
  end if

  call inv_trans(KRESOL=handles(res_idx), PSPSCALAR=spectral_field, PGP=grid_point_field)

  ! -------------- Gather the result on the root (0) and write to file -----------------

  ! Get counts from each PE.
  grid_partition_size_local(1) = num_grid_points
  if (allocated(grid_partition_sizes)) deallocate(grid_partition_sizes)
  allocate(grid_partition_sizes(num_ranks))
  grid_partition_sizes = 0

  call MPI_Gather(grid_partition_size_local, 1, MPI_INT, &
                  grid_partition_sizes, 1, MPI_INT, &
                  0, MPI_COMM_WORLD, ierror)
  if (ierror /= 0) then
    print*,"MPI ERROR => ", ierror
    call ABORT_TRANS("MPI ERROR")
  end if

  print*, "SIZES => ", grid_partition_sizes(:)

  ! Allocate a global field.
  if (allocated(g_grid_point_field)) deallocate(g_grid_point_field)
  allocate(g_grid_point_field(g_num_grid_points))

  ! Make displacement arrays
  if (allocated(displs)) deallocate(displs)
  allocate(displs(num_ranks))
  displs = 0
  do i=2, num_ranks
    displs(i) = displs(i - 1) + grid_partition_sizes(i - 1)
  end do
  print*,"displs => ", displs(:)

  call MPI_Gatherv(grid_point_field(:,1,1), num_grid_points, MPI_FLOAT, &
                   g_grid_point_field, grid_partition_sizes, displs, MPI_FLOAT, &
                   0, MPI_COMM_WORLD, ierror)
  if (ierror /= 0) then
    print*,"Gatherv MPI ERROR => ", ierror
    call ABORT_TRANS("MPI ERROR")
  end if

  if (rank == 1) then
    write(filename, "(A22,I0,A4)") "grid_point_field_trunc", truncation, ".dat"
    open(7, file=filename, form="unformatted")
    write(7) g_grid_point_field(:)
    close(7)
  end if
end do


call MPL_END()
!call MPI_Finalize()

end program
