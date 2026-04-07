/*
 * (C) Crown Copyright 2025- Met Office.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ectrans/transi.h"

#include "transi_test.h"

#include <mpi.h>

// ----------------------------------------------------------------------------

int getColour(const int world_rank) { return world_rank % 2; }

// ----------------------------------------------------------------------------

int main ( int arc, char **argv ) {
  MPI_Init(&arc, &argv);
  trans_use_mpi(true);

  setbuf(stdout,NULL); // unbuffered stdout

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Split world communicator.
  const int colour = getColour(world_rank);
  MPI_Comm split_comm;
  MPI_Comm_split(MPI_COMM_WORLD, colour, world_rank, &split_comm);

  int split_size;
  MPI_Comm_size(split_comm, &split_size);

  // Set default fiat MPL comm.
  const MPI_Fint split_comm_int = MPI_Comm_c2f(split_comm);
  TRANS_CHECK( trans_set_mpi_comm(split_comm_int) );

  // Initialise trans (+ MPL as a result) with split communicator.
  struct Trans_t trans;
  TRANS_CHECK( trans_new(&trans) );

  const int nlon = 320;
  const int nlat = 161;
  const int nsmax = 159;
  TRANS_CHECK( trans_set_resol_lonlat(&trans,nlon,nlat) );
  TRANS_CHECK( trans_set_trunc(&trans,nsmax) );
  TRANS_CHECK( trans_setup(&trans) );

  printf("World size => %d :: Split size => %d :: Trans size => %d\n",
         world_size, split_size, trans.nproc);

  ASSERT_MSG(world_size >= 2,
             "ERROR: Number of MPI processes for this test must be greater than or equal to 2.");
  ASSERT(trans.nproc == split_size);
  ASSERT(trans.nproc < world_size);
  ASSERT(trans.nproc <= world_size / 2);

  // Attempt to set up trans on WORLD - should fail, since MPL has already
  // been initialised on the split_comm
  const MPI_Fint world_comm_int = MPI_Comm_c2f(MPI_COMM_WORLD);
  const int ret_code = trans_set_mpi_comm(world_comm_int);
  ASSERT_MSG(ret_code != 0, "ERROR: Expected `trans_set_mpi_comm(MPI_COMM_WORLD)` "
                            "to fail on second setup.");

  TRANS_CHECK( trans_finalize() );

  MPI_Finalize();

  return 0;
}
