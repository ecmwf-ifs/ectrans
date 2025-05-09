#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ectrans/transi.h"
#include "transi_test.h"

void read_grid(struct Trans_t*);

int main ( int arc, char **argv )
{
  fprintf(stderr,"ectrans version = %s\n",ectrans_version());

  fprintf(stderr,"Using MPI: %d\n", test_use_mpi());

  trans_use_mpi( test_use_mpi() );

  // must be done before trans_new
  trans_set_leq_regions(false);
  //trans_set_nprgpew(1);

  //printf("transi started\n");
  int nout = 3;
  struct Trans_t trans;
  trans_new(&trans);
  read_grid(&trans);
  trans_setup(&trans);

  if( trans.myproc == 1 ) {
    fprintf(stderr,"nproc   = %d\n", trans.nproc);
    fprintf(stderr,"ngptot  = %d\n", trans.ngptot);
    fprintf(stderr,"ngptotg = %d\n", trans.ngptotg);
  }
  //trans_inquire(&trans,"nvalue,mvalue");

  // Allocate gridpoint data
  int nvordiv = 1;
  int nscalar = 2;
  int nfld  = 2*nvordiv+nscalar;
  double* rgp  = malloc( sizeof(double) * nfld *trans.ngptot  );

  // Load data on proc 1
  double* rgpg = NULL;
  if( trans.myproc == 1 )
  {
    rgpg = malloc( sizeof(double) * 4*trans.ngptotg );
    int i, j;
    for ( j=0;j<nvordiv;j++) {
		  for( i=0; i<trans.ngptotg; ++i )
		  {
		    rgpg[(2*j)*trans.ngptotg+i] = 2*j; 		// U
		    rgpg[(2*j+1)*trans.ngptotg+i] = 2*j+1; // V
		  }
		}
    for ( j=0;j<nscalar;j++) {
		  for( i=0; i<trans.ngptotg; ++i )
		  {
		    rgpg[(2*nvordiv+j)*trans.ngptotg+i] =2*nvordiv+j; // scalar
		  }
		}
  }

  // Distribute gridpoint data from proc 1
  int* nfrom = malloc( sizeof(int) * nfld );
  for ( int j=0;j<nfld;j++) nfrom[j]=1;

  struct DistGrid_t distgrid = new_distgrid(&trans);
    distgrid.nfrom = nfrom;
    distgrid.rgpg  = rgpg;
    distgrid.rgp   = rgp;
    distgrid.nfld  = nfld;
  trans_distgrid(&distgrid);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      if(nout < trans.ngptot-1) {
        printf("rgp[%d][...] : ...\n",j);
        i = trans.ngptot-1;
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
    }
  }

  // Allocate spectral data

  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
	double* rmeanu = malloc( sizeof(double) * nvordiv);
	double* rmeanv = malloc( sizeof(double) * nvordiv);
	
  // Direct Transform
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
    dirtrans.nscalar   = nscalar;
    dirtrans.nvordiv   = nvordiv;
    dirtrans.rgp       = rgp;
    dirtrans.rspscalar = rspscalar;
    dirtrans.rspvor    = rspvor;
    dirtrans.rspdiv    = rspdiv;
    dirtrans.rmeanu    = rmeanu;
    dirtrans.rmeanv    = rmeanv;
  trans_dirtrans(&dirtrans);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalar[%d][%d] : %f\n",j,i,rspscalar[i*nscalar+j]);
    }
  }

  // Gather spectral field (for fun)
  int* nto = malloc( sizeof(int) * nscalar );
  nto[0] = 1;
  nto[1] = 1;

  double* rspscalarg = NULL;
  if( trans.myproc == 1 )
    rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );

  struct GathSpec_t gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspscalar;
    gathspec.rspecg = rspscalarg;
    gathspec.nfld   = nscalar;
    gathspec.nto    = nto;
  trans_gathspec(&gathspec);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nvordiv; ++j)
    {
      printf("rmeanu[%d] : %f\n",j,rmeanu[j]);
      printf("rmeanv[%d] : %f\n",j,rmeanv[j]);
      for( i=0; i<nout; ++i )
        printf("rspvor[%d][%d] : %f\n",j,i,rspvor[i*nvordiv+j]);
        printf("rspdiv[%d][%d] : %f\n",j,i,rspdiv[i*nvordiv+j]);
    }
    for( j=0; j<nscalar; ++j)
    {
      for( i=0; i<nout; ++i )
        printf("rspscalarg[%d][%d] : %f\n",j,i,rspscalarg[i*nscalar+j]);
    }
  }
  
  // Distribute spectral field (for fun)
  struct DistSpec_t distspec = new_distspec(&trans);
    distspec.rspec  = rspscalar;
    distspec.rspecg = rspscalarg;
    distspec.nfld   = nscalar;
    distspec.nfrom  = nto;
  trans_distspec(&distspec);

  // Inverse Transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
    invtrans.nscalar   = nscalar;
    invtrans.nvordiv   = nvordiv;
    invtrans.rspscalar = rspscalar;
    invtrans.rspvor    = rspvor;
    invtrans.rspdiv    = rspdiv;
    invtrans.rmeanu    = rmeanu;
    invtrans.rmeanv    = rmeanv;
    invtrans.rgp       = rgp;
  trans_invtrans(&invtrans);

  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
      if(nout < trans.ngptot-1) {
        printf("rgp[%d][...] : ...\n",j);
        i = trans.ngptot-1;
        printf("rgp[%d][%d] : %f\n",j,i,rgp[j*trans.ngptot+i]);
      }
    }
  }

  // Gather gridpoint fields
  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nto  = nfrom;
    gathgrid.nfld = nfld;
  trans_gathgrid(&gathgrid);


  if( trans.myproc == 1 )
  {
    int i,j;
    for( j=0; j<nfld; ++j)
    {
      for( i=0; i<nout; ++i )
      {
        printf("rgpg[%d][%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
      }
      if(nout < trans.ngptot-1) {
        printf("rgpg[%d][...] : ...\n",j);
        i = trans.ngptotg-1;
        printf("rgpg[%d][%d] : %f\n",j,i,rgpg[j*trans.ngptotg+i]);
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
  free(rmeanu);
  free(rmeanv);

  if( trans.myproc == 1 ) {
    printf("cleanup\n");
  }
  trans_delete(&trans);

  trans_finalize();

  // printf("transi finished\n");
  return 0;
}







void read_grid(struct Trans_t* trans)
{
	// lam grid of 20x18
	trans->ndgl=18;
	trans->nlon=20;
	trans->nsmax=9;
	trans->nmsmax=8;
	trans->llam=true;
}
