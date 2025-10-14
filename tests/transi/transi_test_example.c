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

// TODO
#include "/home/users/darth/opt/mobbs-61/util/spice_gnu/include/mpi.h"

static bool check_values = false;

// ----------------------------------------------------------------------------

void test_gptosptogp(const MPI_Comm comm, int nlon, int nlat, int nsmax)
{
  double gptosp_tol = 1.e-6;
  double sptogp_tol = 1.e-6;
  printf("test_gptosptogp( nlon=%d, nlat=%d, nsmax=%d )\n",nlon,nlat,nsmax);

  int nout = 2;

  // Set trans mpi comm
  const MPI_Fint comm_int = MPI_Comm_c2f(comm);
  trans_set_mpi_comm(comm_int);

  struct Trans_t trans;
  TRANS_CHECK( trans_new(&trans) );
  TRANS_CHECK( trans_set_resol_lonlat(&trans,nlon,nlat) );
  TRANS_CHECK( trans_set_trunc(&trans,nsmax) );

  //set_standard_rgg(&trans,(nlat-1)/2,nsmax);


  trans.fft = TRANS_FFTW;
  trans.flt = 0;

  TRANS_CHECK( trans_setup(&trans) );
  printf("ndgl  = %d\n",trans.ndgl);
  printf("nsmax = %d\n",trans.nsmax);
  printf("ngptotg  = %d\n",trans.ngptotg);

  ASSERT(trans.ngptotg == nlon*nlat);

  // In case of odd number of latitudes, there is one latitude extra in distributed gridpoints
  // This can only be checked with nproc==1
  if( trans.nproc == 1 )
    ASSERT(trans.ngptot  == nlon*nlat + nlon*(nlat%2));

  // Allocate gridpoint data
  int nscalar = 2;
  int nvordiv = 1;
  int nfld  = 2*nvordiv+nscalar;
  double* rgp  = malloc( sizeof(double) * nfld * trans.ngptot );

  // Load data on proc 1
  double* rgpg = NULL;
  if( trans.myproc == 1 )
  {
    rgpg = malloc( sizeof(double) * nfld*trans.ngptotg );
    int i;
    for( i=0; i<trans.ngptotg; ++i )
    {
      rgpg[0*trans.ngptotg+i] = 1.; // U
      rgpg[1*trans.ngptotg+i] = 2.; // V
      rgpg[2*trans.ngptotg+i] = 3.; // scalar 1
      rgpg[3*trans.ngptotg+i] = 4.; // scalar 2
    }
  }

  // Distribute gridpoint data from proc 1
  int* nfrom = malloc( sizeof(int) * nfld );
  nfrom[0] = 1; // U
  nfrom[1] = 1; // V
  nfrom[2] = 1; // scalar 1
  nfrom[3] = 1; // scalar 2

  printf("trans_distgrid()\n");
  struct DistGrid_t distgrid = new_distgrid(&trans);
    distgrid.nfrom = nfrom;
    distgrid.rgpg  = rgpg;
    distgrid.rgp   = rgp;
    distgrid.nfld  = nfld;
  TRANS_CHECK( trans_distgrid(&distgrid) );

  if( trans.myproc == 1 )
  {

    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[fld=%d][pt=%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      if ( j>=0 )
      {
        for( i=0; i<trans.ngptot; ++i )
        {
          if( fabs(rgp[j*trans.ngptot+i]-(j+1)) > gptosp_tol)
            printf("error --> rgp[fld=%d][pt=%d] : %f instead of %d\n",j,i,rgp[j*trans.ngptot+i],(j+1));
          ASSERT( fabs(rgp[j*trans.ngptot+i]-(j+1)) < gptosp_tol);
        }
      }

    }
  }

  // Allocate spectral data

  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );

  // Direct Transform
  printf("trans_dirtrans()\n");
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nscalar   = nscalar;
    dirtrans.nvordiv   = nvordiv;
    dirtrans.rgp       = rgp;
    dirtrans.rspscalar = rspscalar;
    dirtrans.rspvor    = rspvor;
    dirtrans.rspdiv    = rspdiv;
  TRANS_CHECK( trans_dirtrans(&dirtrans) );

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalar[fld=%d][wave=%d] : %f\n",j,i,rspscalar[i*nscalar+j]);

      if ( j>=2 && check_values )
      {
        for( i=1; i<trans.nspec2; ++i )
        {
          ASSERT( fabs(rspscalar[i*nscalar+j]) < gptosp_tol );
          if( fabs(rspscalar[i*nscalar+j]) > gptosp_tol )
            printf("error --> rspscalar[fld=%d][wave=%d] : %f\n",j,i,rspscalar[i*nscalar+j]);
        }
      }
    }
  }

  // Gather spectral field (for fun)
  int* nto = malloc( sizeof(int) * nscalar );
  nto[0] = 1;
  nto[1] = 1;

  double* rspscalarg = NULL;
  if( trans.myproc == 1 )
    rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );

  printf("trans_gathspec()\n");
  struct GathSpec_t gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspscalar;
    gathspec.rspecg = rspscalarg;
    gathspec.nfld   = nscalar;
    gathspec.nto    = nto;
  TRANS_CHECK( trans_gathspec(&gathspec) );

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);

      if( check_values )
      {
        for( i=1; i<trans.nspec2g; ++i )
        {
          ASSERT( fabs(rspscalarg[i*nscalar+j]) < gptosp_tol );
          if( fabs(rspscalarg[i*nscalar+j]) > gptosp_tol )
          printf("error -> rspscalarg[fld=%d][wave=%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);
        }
      }
    }
  }

  printf("trans_distspec()\n");
  // Distribute spectral field (for fun)
  struct DistSpec_t distspec = new_distspec(&trans);
    distspec.rspec  = rspscalar;
    distspec.rspecg = rspscalarg;
    distspec.nfld   = nscalar;
    distspec.nfrom  = nto;
  TRANS_CHECK( trans_distspec(&distspec) );


  printf("trans_invtrans()\n");
  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rgp       = rgp;
  TRANS_CHECK( trans_invtrans(&invtrans) );

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      if( j>= 2 && check_values )
      {
        for( i=0; i<trans.ngptot; ++i )
        {
          if ( fabs(rgp[j*trans.ngptot+i] - (j+1) )/(double)(j+1) > sptogp_tol )
            printf("error --> rgp[fld=%d][pt=%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
          ASSERT( fabs(rgp[j*trans.ngptot+i] - (j+1) )/(double)(j+1) < sptogp_tol );
        }
      }
    }
  }

  printf("trans_gathgrid()\n");
  // Gather gridpoint fields
  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nto  = nfrom;
    gathgrid.nfld = nfld;
  TRANS_CHECK( trans_gathgrid(&gathgrid) );


  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgpg[fld=%d][pt=%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
      }
      if( j>=2 && check_values )
      {
        for( i=0; i<trans.ngptotg; ++i )
        {
          ASSERT( fabs(rgpg[j*trans.ngptotg+i] - (j+1) )/(double)(j+1) < sptogp_tol );
          if ( fabs(rgpg[j*trans.ngptotg+i] - (j+1) )/(double)(j+1) > sptogp_tol )
            printf("error --> rgpg[fld=%d][pt=%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
        }
      }
    }
  }

  // Deallocate arrays
  free(rgp);
  free(rgpg);
  free(rspscalar);
  free(rspscalarg);
  free(rspvor);
  free(rspdiv);
  free(nfrom);
  free(nto);

  TRANS_CHECK( trans_delete(&trans) );
}

// ----------------------------------------------------------------------------

int main ( int arc, char **argv )
{
  MPI_Init(&arc, &argv);
  trans_use_mpi(true);

  setbuf(stdout,NULL); // unbuffered stdout

  printf("-----------------------------\n");

  // Split for IO
  MPI_Comm base_comm;
  bool is_io_server = false;
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    is_io_server = (rank == 0);
    const int colour = (int)is_io_server;

    MPI_Comm_split(MPI_COMM_WORLD, colour, rank, &base_comm);
  }

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (!is_io_server) {
    // Split the comm
    int colour = -1;
    MPI_Comm split_comm;
    {
      int world_size, rank;
      MPI_Comm_size(base_comm, &world_size);
      MPI_Comm_rank(base_comm, &rank);

      const float ratio = 0.25;  // 1:3 rank splitting.
      if ((float)rank/(float)world_size >= ratio) {
        colour = 0;
      } else {
        colour = 1;
      }

      MPI_Comm_split(base_comm, colour, rank, &split_comm);
    }

    if (colour == 0) {
      printf("WORLD_RANK %d, g%d => L\n", world_rank, colour);
      test_gptosptogp(split_comm, 320,161,159);
    } else if (colour == 1) {
      printf("WORLD_RANK %d, g%d => S\n", world_rank, colour);
      test_gptosptogp(split_comm, 320/2,161/2,159/2);
    }

    TRANS_CHECK( trans_finalize() );
  }

  printf("WORLD_RANK %d => ~Done~\n", world_rank);

  MPI_Finalize();

  return 0;
}


