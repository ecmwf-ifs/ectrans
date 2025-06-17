#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ectrans/transi.h"
#include "transi_test.h"

#define print(...) fprintf(stderr, __VA_ARGS__)

void print_gp(const char* name, double* array, int nfld, int ngp, int nout);
void print_sp(const char* name, double* array, int nfld, int nspec2, int nout);

int main ( int arc, char **argv ) {

  print("ectrans version = %s\n",ectrans_version());
  print("Using MPI: %d\n", test_use_mpi());

  // Must be done before first trans_setup or trans_init
  trans_use_mpi( test_use_mpi() );
  trans_set_leq_regions(false);
  trans_set_nprgpew(transi_test_nprgpew());

  int nout_gp = 3;
  int nout_sp = 3;
  struct Trans_t trans;
  trans_new(&trans);
  // lam grid of 20x18
  int nx = 20;
  int ny = 18;
  double dx = 2500.0;
  double dy = 2500.0;
  int tx = (nx-1)/2;
  int ty = (ny-1)/2;
  trans_set_resol_lam(&trans, nx, ny, dx, dy);
  trans_set_trunc_lam(&trans, ty, tx);
  trans_setup(&trans);

  if( trans.myproc == 1 ) {
    print("nproc   = %d\n", trans.nproc);
    print("nx , ny = %d , %d \n", nx, ny);
    print("tx , ty = %d , %d \n", tx, ty);
    print("ngptot  = %d\n", trans.ngptot);
    print("ngptotg = %d\n", trans.ngptotg);
  }

  // Allocate gridpoint data
  int nvordiv = 1;
  int nscalar = 2;
  int nfld  = 2*nvordiv+nscalar;
  double* rgp  = malloc( sizeof(double) * nfld *trans.ngptot  );

  // Load data on proc 1
  double* rgpg = NULL;
  if( trans.myproc == 1 ) {
    rgpg = malloc( sizeof(double) * 4*trans.ngptotg );
    for (int j=0;j<nvordiv;j++) {
      for(int i=0; i<trans.ngptotg; ++i) {
        rgpg[(2*j+0)*trans.ngptotg+i] = 2*j+1; // U
        rgpg[(2*j+1)*trans.ngptotg+i] = 2*j+2; // V
      }
    }
    for (int j=0;j<nscalar;j++) {
      for(int i=0; i<trans.ngptotg; ++i) {
        rgpg[(2*nvordiv+j)*trans.ngptotg+i] = 2*nvordiv+j + 1; // scalar
      }
    }
  }

  // Distribute gridpoint data from proc 1
  int* nfrom = malloc( sizeof(int) * nfld );
  for (int j=0;j<nfld;j++) nfrom[j]=1;

  struct DistGrid_t distgrid = new_distgrid(&trans);
    distgrid.nfrom = nfrom;
    distgrid.rgpg  = rgpg;
    distgrid.rgp   = rgp;
    distgrid.nfld  = nfld;
  trans_distgrid(&distgrid);

  if( trans.myproc == 1 ) {
    print_gp("rgp", rgp, nfld, trans.ngptot, nout_gp);
  }

  // Allocate spectral data

  double* rspvor    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspdiv    = malloc( sizeof(double) * nvordiv*trans.nspec2 );
  double* rspscalar = malloc( sizeof(double) * nscalar*trans.nspec2 );
  double* rmeanu    = malloc( sizeof(double) * nvordiv);
  double* rmeanv    = malloc( sizeof(double) * nvordiv);

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

  if( trans.myproc == 1 ) {
    for(int j=0; j<nvordiv; ++j) {
      print("rmeanu[%d] : %f\n",j,rmeanu[j]);
      print("rmeanv[%d] : %f\n",j,rmeanv[j]);
    }
    print_sp("rspvor",rspvor,nvordiv,trans.nspec,nout_sp);
    print_sp("rspdiv",rspdiv,nvordiv,trans.nspec,nout_sp);
    print_sp("rspscalar",rspscalar,nscalar,trans.nspec,nout_sp);
  }

  // Gather spectral field (for fun)
  int* nto = malloc( sizeof(int) * nscalar );
  nto[0] = 1;
  nto[1] = 1;

  double* rspscalarg = NULL;
  if( trans.myproc == 1 ) {
    rspscalarg = malloc( sizeof(double) * nscalar*trans.nspec2g );
  }

  struct GathSpec_t gathspec = new_gathspec(&trans);
    gathspec.rspec  = rspscalar;
    gathspec.rspecg = rspscalarg;
    gathspec.nfld   = nscalar;
    gathspec.nto    = nto;
  trans_gathspec(&gathspec);

  if( trans.myproc == 1 ) {
    print_sp("rspscalarg",rspscalarg,nscalar,trans.nspec2g,nout_sp);
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

  if( trans.myproc == 1 ) {
    print_gp("rgp",rgp,nfld,trans.ngptot,nout_gp);
  }

  // Gather gridpoint fields
  struct GathGrid_t gathgrid = new_gathgrid(&trans);
    gathgrid.rgp  = rgp;
    gathgrid.rgpg = rgpg;
    gathgrid.nto  = nfrom;
    gathgrid.nfld = nfld;
  trans_gathgrid(&gathgrid);


  if( trans.myproc == 1 ) {
    print_gp("rgpg",rgpg,nfld,trans.ngptotg,nout_gp);
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
    print("cleanup\n");
  }
  trans_delete(&trans);

  trans_finalize();

  return 0;
}

// ---------------------------------------------------------------------------------

void print_gp(const char* name, double* array, int nfld, int ngp, int nout) {
  if (nout < 0) {
    nout = ngp;
  }
  for(int j=0; j<nfld; ++j) {
    for(int i=0; i<nout; ++i) {
      print("%s[%d][%3d] : %f\n",name,j,i,array[j*ngp+i]);
    }
    if(nout < ngp-1) {
      print("%s[%d][...] : ...\n",name,j);
      int i = ngp-1;
      print("%s[%d][%3d] : %f\n",name,j,i,array[j*ngp+i]);
    }
  }
}

void print_sp(const char* name, double* array, int nfld, int nspec2, int nout) {
  if (nout < 0) {
    nout = nspec2;
  }
  for(int j=0; j<nfld; ++j) {
    for(int i=0; i<nout; ++i ) {
      print("%s[%3d][%d] : %f\n",name, i,j,array[i*nfld+j]);
    }
    if(nout < nspec2-1) {
      print("%s[...][%d] : ...\n",name,j);
      int i = nspec2-1;
      print("%s[%3d][%d] : %f\n",name, i,j,array[i*nfld+j]);
    }
  }
}
