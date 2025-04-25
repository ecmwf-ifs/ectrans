/*
 * (C) Crown Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ectrans/transi.h"
#include "transi_test.h"

// ----------------------------------------------------------------------------

void test_invtrans_adjoint(int nlon, int nlat, int nsmax)
{
const unsigned int seed = 123;
double adjoint_tol = 1.e-12;
printf("test_invtrans_adjoint( nlon=%d, nlat=%d, nsmax=%d )\n",nlon,nlat,nsmax);

// ===== Set-up trans =====
struct Trans_t trans;
TRANS_CHECK( trans_new(&trans) );

int* nloen  = malloc( sizeof(int) * nlat);
{
    int i;
    for( i=0; i<nlat; ++i )
      nloen[i] = nlon;
}
TRANS_CHECK( trans_set_resol(&trans,nlat, nloen) );
TRANS_CHECK( trans_set_trunc(&trans,nsmax) );

trans.fft = TRANS_FFTW;
trans.flt = 0;

TRANS_CHECK( trans_setup(&trans) );
printf("ndgl  = %d\n",trans.ndgl);
printf("nsmax = %d\n",trans.nsmax);
printf("ngptotg  = %d\n",trans.ngptotg);
printf("llatlon  = %d\n",trans.llatlon);
ASSERT(trans.ngptotg == nlon*nlat);

int nscalar = 2;
int nvordiv = 1;
int nfld  = 2*nvordiv+nscalar;

// ===== Allocate and initialize gridpoint data =====
double* rgpx  = calloc( nfld * trans.ngptot, sizeof(double) );
double* rgpy  = calloc( nfld * trans.ngptot, sizeof(double) );
double* rgpxg = NULL;
double* rgpyg = NULL;
if( trans.myproc == 1 ) {
  rgpxg = malloc( sizeof(double) * nfld * trans.ngptotg );
  rgpyg = calloc( nfld * trans.ngptotg, sizeof(double) );
}

// Initialize global spectral data on proc 1
srand(seed);
if( trans.myproc == 1 )
{
  int i;
  for( i=0; i<trans.ngptotg; ++i )
  {
    rgpxg[0*trans.ngptotg+i] = rand() * (2.0 / RAND_MAX) - 1.0; // U
    rgpxg[1*trans.ngptotg+i] = rand() * (2.0 / RAND_MAX) - 1.0; // V
    rgpxg[2*trans.ngptotg+i] = rand() * (2.0 / RAND_MAX) - 1.0; // scalar 1
    rgpxg[3*trans.ngptotg+i] = rand() * (2.0 / RAND_MAX) - 1.0; // scalar 2
  }
}

// Distribute gridpoint data from proc 1 to initialize local gridpoint data
int* nfrom = malloc( sizeof(int) * nfld );
nfrom[0] = 1; // U
nfrom[1] = 1; // V
nfrom[2] = 1; // scalar 1
nfrom[3] = 1; // scalar 2

printf("trans_distgrid()\n");
struct DistGrid_t distgrid = new_distgrid(&trans);
  distgrid.nfrom = nfrom;
  distgrid.rgpg  = rgpxg;
  distgrid.rgp   = rgpx;
  distgrid.nfld  = nfld;
TRANS_CHECK( trans_distgrid(&distgrid) );

// ===== Allocate and initialize spectral data =====
double* rspscalarx = calloc( nscalar*trans.nspec2, sizeof(double));
double* rspscalary = calloc( nscalar*trans.nspec2, sizeof(double));
double* rspvorx    = calloc( nvordiv*trans.nspec2, sizeof(double));
double* rspvory    = calloc( nvordiv*trans.nspec2, sizeof(double));
double* rspdivx    = calloc( nvordiv*trans.nspec2, sizeof(double));
double* rspdivy    = calloc( nvordiv*trans.nspec2, sizeof(double));
double* rspscalarxg = NULL;
double* rspscalaryg = NULL;
double* rspvorxg = NULL;
double* rspvoryg = NULL;
double* rspdivxg = NULL;
double* rspdivyg = NULL;
if( trans.myproc == 1 ) {
  rspscalarxg = malloc( sizeof(double) * nscalar * trans.nspec2g );
  rspscalaryg = calloc( nscalar * trans.nspec2g, sizeof(double) );
  rspvorxg = malloc( sizeof(double) * nvordiv * trans.nspec2g );
  rspvoryg = calloc( nvordiv * trans.nspec2g, sizeof(double) );
  rspdivxg = malloc( sizeof(double) * nvordiv * trans.nspec2g );
  rspdivyg = calloc( nvordiv * trans.nspec2g, sizeof(double) );
}

// Initialize global spectral data on proc 1
if( trans.myproc == 1 )
{
  int i;
  for( i=0; i<trans.nspec2g; ++i )
  {
    rspvorxg[i] = rand() * (2.0 / RAND_MAX) - 1.0; // vor
    rspdivxg[i] = rand() * (2.0 / RAND_MAX) - 1.0; // div
    rspscalarxg[0*trans.nspec2g+i] = rand() * (2.0 / RAND_MAX) - 1.0; // scalar 1
    rspscalarxg[1*trans.nspec2g+i] = rand() * (2.0 / RAND_MAX) - 1.0; // scalar 2
  }
}

// Distribute global spectral data from proc 1 to initialize local spectral data
int* nto = malloc( sizeof(int) * nscalar );
nto[0] = 1; // scalar 1 / vor
nto[1] = 1; // scalar 2 / div

printf("trans_distspec()\n");
struct DistSpec_t distspec = new_distspec(&trans);
  distspec.rspec  = rspscalarx;
  distspec.rspecg = rspscalarxg;
  distspec.nfld   = nscalar;
  distspec.nfrom  = nto;
TRANS_CHECK( trans_distspec(&distspec) );
distspec = new_distspec(&trans);
  distspec.rspec  = rspvorx;
  distspec.rspecg = rspvorxg;
  distspec.nfld   = nvordiv;
  distspec.nfrom  = nto;
TRANS_CHECK( trans_distspec(&distspec) );
distspec = new_distspec(&trans);
  distspec.rspec  = rspdivx;
  distspec.rspecg = rspdivxg;
  distspec.nfld   = nvordiv;
  distspec.nfrom  = nto;
TRANS_CHECK( trans_distspec(&distspec) );

// ===== Compute adjoint invtrans and gather result on proc 1 =====
// i.e. invtrans_adj(rgpx) = (rspscalary, rspvory, rspdivy)
printf("trans_invtrans_adj()\n");
struct InvTransAdj_t invtrans_adj = new_invtrans_adj(&trans);
  invtrans_adj.nscalar   = nscalar;
  invtrans_adj.nvordiv   = nvordiv;
  invtrans_adj.rgp       = rgpx;
  invtrans_adj.rspscalar = rspscalary;
  invtrans_adj.rspvor    = rspvory;
  invtrans_adj.rspdiv    = rspdivy;
TRANS_CHECK( trans_invtrans_adj(&invtrans_adj) );

printf("trans_gathspec()\n");
struct GathSpec_t gathspec = new_gathspec(&trans);
  gathspec.rspec  = rspscalary;
  gathspec.rspecg = rspscalaryg;
  gathspec.nfld   = nscalar;
  gathspec.nto    = nto;
TRANS_CHECK( trans_gathspec(&gathspec) );
gathspec = new_gathspec(&trans);
  gathspec.rspec  = rspvory;
  gathspec.rspecg = rspvoryg;
  gathspec.nfld   = nvordiv;
  gathspec.nto    = nto;
TRANS_CHECK( trans_gathspec(&gathspec) );
gathspec = new_gathspec(&trans);
  gathspec.rspec  = rspdivy;
  gathspec.rspecg = rspdivyg;
  gathspec.nfld   = nvordiv;
  gathspec.nto    = nto;
TRANS_CHECK( trans_gathspec(&gathspec) );

// ===== Compute: adj_value1 = <invtrans_adj(rgpx), (rspscalarx, rspvorx, rspdivx)> =====
// i.e. adj_value1 = <(rspscalary, rspvory, rspdivy), (rspscalarx, rspvorx, rspdivx)>

double adj_value1 = 0.0;
if( trans.myproc == 1 )
{
  int i,j;

  for( j=0; j<nscalar; ++j)
  {
    // The first 2*nlat terms are on m=0
    for( i=0; i<2*nlat; ++i )
      adj_value1 += rspscalaryg[i*nscalar+j] * rspscalarxg[i*nscalar+j];

    for( i=2*nlat; i<trans.nspec2g; ++i )
      adj_value1 += 2.0 * rspscalaryg[i*nscalar+j] * rspscalarxg[i*nscalar+j];

    for( i=0; i<trans.nspec2g; i+=2 ) {
      printf("rspscalaryg[fld=%d][coeff=%d].real : %f\n",j,i/2,rspscalaryg[i*nscalar+j]);
      printf("rspscalaryg[fld=%d][coeff=%d].imag : %f\n",j,i/2,rspscalaryg[(i+1)*nscalar+j]);
    }
  }
  
  for( j=0; j<nvordiv; ++j)
  {
    // The first 2*nlat terms are on m=0
    for( i=0; i<2*nlat; ++i ) {
      adj_value1 += rspvoryg[i*nvordiv+j] * rspvorxg[i*nvordiv+j];
      adj_value1 += rspdivyg[i*nvordiv+j] * rspdivxg[i*nvordiv+j];
    }
    
    for( i=2*nlat; i<trans.nspec2g; ++i ) {
      adj_value1 += 2.0 * rspvoryg[i*nvordiv+j] * rspvorxg[i*nvordiv+j];
      adj_value1 += 2.0 * rspdivyg[i*nvordiv+j] * rspdivxg[i*nvordiv+j];
    }

    for( i=0; i<trans.nspec2g; i+=2 ) {
      printf("rspvoryg[fld=%d][coeff=%d].real : %f\n",j,i/2,rspvoryg[i*nvordiv+j]);
      printf("rspvoryg[fld=%d][coeff=%d].imag : %f\n",j,i/2,rspvoryg[(i+1)*nvordiv+j]);
    }
    for( i=0; i<trans.nspec2g; i+=2 ) {
      printf("rspdivyg[fld=%d][coeff=%d].real : %f\n",j,i/2,rspdivyg[i*nvordiv+j]);
      printf("rspdivyg[fld=%d][coeff=%d].imag : %f\n",j,i/2,rspdivyg[(i+1)*nvordiv+j]);
    }
  }
}

// ===== Compute invtrans and gather result on proc 1 =====
// i.e. invtrans(rspscalarx, rspvorx, rspdivx) = rgpy
printf("trans_invtrans()\n");
struct InvTrans_t invtrans = new_invtrans(&trans);
  invtrans.nscalar   = nscalar;
  invtrans.nvordiv   = nvordiv;
  invtrans.rspscalar = rspscalarx;
  invtrans.rspvor    = rspvorx;
  invtrans.rspdiv    = rspdivx;
  invtrans.rgp       = rgpy;
TRANS_CHECK( trans_invtrans(&invtrans) );

printf("trans_gathgrid()\n");
struct GathGrid_t gathgrid = new_gathgrid(&trans);
gathgrid.rgp  = rgpy;
gathgrid.rgpg = rgpyg;
gathgrid.nto  = nfrom;
gathgrid.nfld = nfld;
TRANS_CHECK( trans_gathgrid(&gathgrid) );


// ===== Compute: adj_value2 = <rgpx, invtrans(rspscalarx, rspvorx, rspdivx)> =====
// i.e. adj_value2 = <rgpx, rgpy>
double adj_value2 = 0.0;
if( trans.myproc == 1 )
{
  int i,j;
  for( j=0; j<nfld; ++j)
  { 
    if( j>=0 )
    {
      for( i=0; i<trans.ngptotg; ++i )
        adj_value2 += rgpxg[j*trans.ngptotg+i] * rgpyg[j*trans.ngptotg+i];
    }

    for( i=0; i<trans.ngptotg; ++i ) {
      printf("rgpyg[fld=%d][pt=%d] : %f\n",j,i,rgpyg[j*trans.ngptotg+i]);
    }
  }
}

// ===== Compare inner products =====
// i.e. <invtrans_adj(rgpx), (rspscalarx, rspvorx, rspdivx)> == <rgpx, invtrans(rspscalarx, rspvorx, rspdivx)>
if( trans.myproc == 1 ) {
  printf("[adjval1=%f][adjval2=%f] :\n", adj_value1, adj_value2);
  ASSERT( fabs(adj_value1 - adj_value2 )/adj_value1 < adjoint_tol );
}

// ===== Deallocate arrays and clean up trans =====
free(nloen);
free(rgpx);
free(rgpy);
free(rgpxg);
free(rgpyg);
free(rspscalarx);
free(rspscalary);
free(rspscalarxg);
free(rspscalaryg);
free(rspvorx);
free(rspvory);
free(rspvorxg);
free(rspvoryg);
free(rspdivx);
free(rspdivy);
free(rspdivxg);
free(rspdivyg);
free(nfrom);
free(nto);

TRANS_CHECK( trans_delete(&trans) );

}
//
// ----------------------------------------------------------------------------

int main ( int arc, char **argv )
{
  trans_use_mpi( test_use_mpi() );

  setbuf(stdout,NULL); // unbuffered stdout

// The adjoint test works for standard gaussian latitude grid
//   with no points on the equator or poles.
// nsmax = nlat - 1

  printf("-----------------------------\n");
  test_invtrans_adjoint(8,4,3);

  TRANS_CHECK( trans_finalize() );

  return 0;
}

