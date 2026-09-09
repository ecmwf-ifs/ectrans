/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>

#include "ectrans/transi.h"
#include "transi_test.h"

int main( int argc, char** argv )
{
  (void) argc;
  (void) argv;

  setbuf( stderr, NULL );

  trans_use_mpi( test_use_mpi() );

  struct Trans_t trans;
  TRANS_CHECK( trans_new( &trans ) );
  TRANS_CHECK( trans_set_resol_lonlat( &trans, 16, 9 ) );
  TRANS_CHECK( trans_set_trunc( &trans, 2 ) );
  TRANS_CHECK( trans_setup( &trans ) );

  enum { nfld = 2 };
  int has_zero_local_spectral = trans.nspec2 == 0;
  double* rspec = has_zero_local_spectral ? NULL : calloc( nfld * trans.nspec2, sizeof( double ) );
  double* rspecg = NULL;
  double rnorm[nfld] = {0., 0.};
  int owner[nfld] = {1, 1};

  if( trans.myproc == 1 )
    rspecg = calloc( nfld * trans.nspec2g, sizeof( double ) );

  if( !has_zero_local_spectral && rspec == NULL ) {
    fprintf( stderr, "rank %d/%d: rspec allocation failed\n", trans.myproc, trans.nproc );
    return 1;
  }

  if( trans.myproc == 1 && rspecg == NULL ) {
    fprintf( stderr, "rank %d/%d: rspecg allocation failed\n", trans.myproc, trans.nproc );
    return 1;
  }

  fprintf( stderr,
           "rank %d/%d: nspec=%d nspec2=%d nspec2g=%d nspec2mx=%d nump=%d nfld=%d rspec=%p rspecg=%p\n",
           trans.myproc, trans.nproc, trans.nspec, trans.nspec2, trans.nspec2g, trans.nspec2mx,
           trans.nump, nfld, (void*) rspec, (void*) rspecg );

  if( trans.myproc == trans.nproc )
    ASSERT( has_zero_local_spectral );

  fprintf( stderr, "rank %d/%d: trans_gathspec begin\n", trans.myproc, trans.nproc );
  struct GathSpec_t gathspec = new_gathspec( &trans );
    gathspec.rspec = rspec;
    gathspec.rspecg = rspecg;
    gathspec.nfld = nfld;
    gathspec.nto = owner;
  TRANS_CHECK( trans_gathspec( &gathspec ) );
  fprintf( stderr, "rank %d/%d: trans_gathspec end\n", trans.myproc, trans.nproc );

  fprintf( stderr, "rank %d/%d: trans_distspec begin\n", trans.myproc, trans.nproc );
  struct DistSpec_t distspec = new_distspec( &trans );
    distspec.rspec = rspec;
    distspec.rspecg = rspecg;
    distspec.nfld = nfld;
    distspec.nfrom = owner;
  TRANS_CHECK( trans_distspec( &distspec ) );
  fprintf( stderr, "rank %d/%d: trans_distspec end\n", trans.myproc, trans.nproc );

  fprintf( stderr, "rank %d/%d: trans_specnorm begin\n", trans.myproc, trans.nproc );
  struct SpecNorm_t specnorm = new_specnorm( &trans );
    specnorm.rspec = rspec;
    specnorm.rnorm = rnorm;
    specnorm.nfld = nfld;
    specnorm.nmaster = 1;
  TRANS_CHECK( trans_specnorm( &specnorm ) );
  fprintf( stderr, "rank %d/%d: trans_specnorm end\n", trans.myproc, trans.nproc );

  free( rspec );
  free( rspecg );
  TRANS_CHECK( trans_delete( &trans ) );
  TRANS_CHECK( trans_finalize() );

  return 0;
}
