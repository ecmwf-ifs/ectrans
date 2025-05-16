#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "ectrans/transi.h"
#include "transi_test.h"

// ----------------------------------------------------------------------------

double randomDouble() {
    uint64_t r53 = ((uint64_t)(rand()) << 21) ^ (rand() >> 2);
    return (double)r53 / 9007199254740991.0; // 2^53 - 1
}

// ----------------------------------------------------------------------------

void test_lam_invtrans_adjoint() {

double adjoint_tol = 1.e-6;
printf("test_lam_invtrans_adjoint()\n");

struct Trans_t trans;
TRANS_CHECK( trans_new(&trans) );

TRANS_CHECK( trans_set_resol_lam(&trans, 20, 18) );
TRANS_CHECK( trans_set_trunc_lam(&trans, 9, 8) );
TRANS_CHECK( trans_setup(&trans) );
TRANS_CHECK( trans_inquire(&trans,"nvalue,mvalue") );

// Number of fields
int nscalar = 2;
int nvordiv = 1;
int nfld  = 2*nvordiv+nscalar;

// Allocate test data
double* rgp  = malloc( sizeof(double) * nfld * trans.ngptot );
double* rspvor = malloc( sizeof(double) * nvordiv * trans.nspec2 );
double* rspdiv = malloc( sizeof(double) * nvordiv * trans.nspec2 );
double* rspscalar = malloc( sizeof(double) * nscalar * trans.nspec2 );
double* rmeanu = calloc( nvordiv, sizeof(double));
double* rmeanv = calloc( nvordiv, sizeof(double));
double* rgp1  = malloc( sizeof(double) * nfld * trans.ngptot );
double* rspvor1 = malloc( sizeof(double) * nvordiv * trans.nspec2 );
double* rspdiv1 = malloc( sizeof(double) * nvordiv * trans.nspec2 );
double* rspscalar1 = malloc( sizeof(double) * nscalar * trans.nspec2 );
double* rmeanu1 = calloc( nvordiv, sizeof(double));
double* rmeanv1 = calloc( nvordiv, sizeof(double));
double* rgp2 = calloc( nfld * trans.ngptot, sizeof(double));
double* rspvor2 = calloc( nvordiv * trans.nspec2, sizeof(double));
double* rspdiv2 = calloc( nvordiv * trans.nspec2, sizeof(double));
double* rspscalar2 = calloc( nscalar * trans.nspec2, sizeof(double));
double* rmeanu2 = calloc( nvordiv, sizeof(double));
double* rmeanv2 = calloc( nvordiv, sizeof(double));
int* nfrom = malloc( sizeof(int) * nfld );

// Data distribution vector
nfrom[0] = 1; // U
nfrom[1] = 1; // V
nfrom[2] = 1; // scalar 1
nfrom[3] = 1; // scalar 2

// Create random grid-point data on proc 1
double* rgpg = NULL;
if ( trans.myproc == 1 ) {
  rgpg = malloc( sizeof(double) * nfld * trans.ngptotg );
  for( int i=0; i<trans.ngptotg; ++i ) {
    rgpg[0*trans.ngptotg+i] = randomDouble(); // U
    rgpg[1*trans.ngptotg+i] = randomDouble(); // V
    rgpg[2*trans.ngptotg+i] = randomDouble(); // scalar 1
    rgpg[3*trans.ngptotg+i] = randomDouble(); // scalar 2
  }
}

// Distribute gridpoint data from proc 1
struct DistGrid_t distgrid = new_distgrid(&trans);
  distgrid.nfrom = nfrom;
  distgrid.rgpg  = rgpg;
  distgrid.rgp   = rgp1;
  distgrid.nfld  = nfld;
TRANS_CHECK( trans_distgrid(&distgrid) );

// Create random spectral data on proc 1
double* rspg = NULL;
if ( trans.myproc == 1 ) {
  rspg = calloc( nfld * trans.nspec2g, sizeof(double) );
  for( int i=0; i<trans.nspec2g; ++i ) {
    rspg[i*nfld+0] = randomDouble(); // vorticity
    rspg[i*nfld+1] = randomDouble(); // divergence
    rspg[i*nfld+2] = randomDouble(); // scalar 1
    rspg[i*nfld+3] = randomDouble(); // scalar 2
  }
}

// Distribute spectral data from proc 1
double* rsp  = malloc( sizeof(double) * nfld * trans.nspec2 );
struct DistSpec_t distspec = new_distspec(&trans);
  distspec.rspec  = rsp;
  distspec.rspecg = rspg;
  distspec.nfld   = nfld;
  distspec.nfrom  = nfrom;
TRANS_CHECK( trans_distspec(&distspec) );

// Split spectral fields
for( int i=0; i<trans.nspec2; ++i ) {
  rspvor2[i*nvordiv+0] = rsp[i*nfld+0]; // vorticity
  rspdiv2[i*nvordiv+0] = rsp[i*nfld+1]; // divergence
  rspscalar2[i*nscalar+0] = rsp[i*nfld+2]; // scalar 1
  rspscalar2[i*nscalar+1] = rsp[i*nfld+3]; // scalar 2
}

// Initialize mean wind
for( int j=0; j<nvordiv; ++j) {
  rmeanu2[j] = randomDouble();
  rmeanv2[j] = randomDouble();
}

// Apply inverse transform
struct InvTrans_t invtrans = new_invtrans(&trans);
  invtrans.nscalar   = nscalar;
  invtrans.nvordiv   = nvordiv;
  invtrans.rspscalar = rspscalar2;
  invtrans.rspvor    = rspvor2;
  invtrans.rspdiv    = rspdiv2;
  invtrans.rmeanu    = rmeanu2;
  invtrans.rmeanv    = rmeanv2;
  invtrans.rgp       = rgp2;
TRANS_CHECK( trans_invtrans(&invtrans) );

// Compute grid-point dot product
double adj_value_gp = 0.0;
for(int j=0; j<nfld; ++j)
  for( int i=0; i<trans.ngptot; ++i )
      adj_value_gp += rgp1[j*trans.ngptot+i] * rgp2[j*trans.ngptot+i];
printf("Grid-point product = %.12f\n", adj_value_gp);

// Apply adjoint of the direct transform
struct InvTransAdj_t invtrans_adj = new_invtrans_adj(&trans);
  invtrans_adj.nscalar   = nscalar;
  invtrans_adj.nvordiv   = nvordiv;
  invtrans_adj.rgp       = rgp1;
  invtrans_adj.rspscalar = rspscalar1;
  invtrans_adj.rspvor    = rspvor1;
  invtrans_adj.rspdiv    = rspdiv1;
  invtrans_adj.rmeanu    = rmeanu1;
  invtrans_adj.rmeanv    = rmeanv1;
TRANS_CHECK( trans_invtrans_adj(&invtrans_adj) );

// Compute spectral dot product
double adj_value_sp = 0.0;
for( int i=0; i<trans.nspec2; ++i ) {
  adj_value_sp += rspvor1[i] * rspvor2[i];
  adj_value_sp += rspdiv1[i] * rspdiv2[i];
}
for( int j=0; j<nscalar; ++j)
  for( int i=0; i<trans.nspec2; ++i ) {
    adj_value_sp += rspscalar1[i*nscalar+j] * rspscalar2[i*nscalar+j];
}
for( int j=0; j<nvordiv; ++j) {
  adj_value_sp += rmeanu1[j] * rmeanu2[j];
  adj_value_sp += rmeanv1[j] * rmeanv2[j];
}
printf("Spectral product   = %.12f\n", adj_value_sp);

// Compare dot products
ASSERT( fabs(adj_value_sp - adj_value_gp ) / fabs( adj_value_gp ) < adjoint_tol );

// Deallocate arrays
free(nfrom);
free(rgpg);
free(rgp1);
free(rspg);
free(rsp);
free(rspscalar1);
free(rspvor1);
free(rspdiv1);
free(rgp2);
free(rspscalar2);
free(rspvor2);
free(rspdiv2);

TRANS_CHECK( trans_delete(&trans) );

}

// ----------------------------------------------------------------------------

int main ( int arc, char **argv ) {
  trans_use_mpi( test_use_mpi() );

  setbuf(stdout,NULL); // unbuffered stdout

// The adjoint test works for standard gaussian latitude grid
//   with no points on the equator or poles.
// nsmax = nlat - 1

  printf("-----------------------------\n");
  test_lam_invtrans_adjoint();

  TRANS_CHECK( trans_finalize() );

  return 0;
}
