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
#include <math.h>
#include <string.h>

#include "ectrans/transi.h"
#include "transi_test.h"

static int spectral_index( int coeff, int field, int nfld ) {
  return coeff * nfld + field;
}

static double spectral_value( int field, int m, int n, int part ) {
  return 10. * ( field + 1 ) + 3. * ( m + 1 ) + 0.25 * ( n + 1 ) + 0.125 * part;
}

static void setup_global_spectrum( double* rspecg, int nfld, int nsmax ) {
  int coeff = 0;
  for( int m = 0; m <= nsmax; ++m ) {
    for( int n = m; n <= nsmax; ++n ) {
      for( int field = 0; field < nfld; ++field ) {
        rspecg[spectral_index( coeff, field, nfld )] = spectral_value( field, m, n, 0 );
        rspecg[spectral_index( coeff + 1, field, nfld )] = m == 0 ? 0. : spectral_value( field, m, n, 1 );
      }
      coeff += 2;
    }
  }
}

static void compute_global_spectral_norm_from_spectral_values( int nfld, int nsmax, double* rnorm ) {
  for( int field = 0; field < nfld; ++field ) {
    rnorm[field] = 0.;
  }

  for( int m = 0; m <= nsmax; ++m ) {
    for( int n = m; n <= nsmax; ++n ) {
      for( int field = 0; field < nfld; ++field ) {
        const double real = spectral_value( field, m, n, 0 );
        const double imag = m == 0 ? 0. : spectral_value( field, m, n, 1 );
        rnorm[field] += m == 0 ? real * real : 2. * ( real * real + imag * imag );
      }
    }
  }

  for( int field = 0; field < nfld; ++field ) {
    rnorm[field] = sqrt( rnorm[field] );
  }
}

static void compute_global_spectral_norm_from_rspecg( const double* rspecg, int nfld, int nsmax, double* rnorm ) {
  for( int field = 0; field < nfld; ++field ) {
    rnorm[field] = 0.;
  }

  int coeff = 0;
  for( int m = 0; m <= nsmax; ++m ) {
    for( int n = m; n <= nsmax; ++n ) {
      (void) n; // avoid warning of unused variable
      for( int field = 0; field < nfld; ++field ) {
        const double real = rspecg[spectral_index( coeff, field, nfld )];
        const double imag = rspecg[spectral_index( coeff + 1, field, nfld )];
        rnorm[field] += m == 0 ? real * real : 2. * ( real * real + imag * imag );
      }
      coeff += 2;
    }
  }

  for( int field = 0; field < nfld; ++field )
    rnorm[field] = sqrt( rnorm[field] );
}

static int is_approx_eq( double actual, double expected ) {
  const double tolerance = 1.e-12 * fmax( 1., fabs( expected ) );
  return fabs( actual - expected ) <= tolerance;
}

static int spectra_are_approx_eq( const double* actual, const double* expected, int size ) {
  for( int i = 0; i < size; ++i ) {
    if( !is_approx_eq( actual[i], expected[i] ) )
      return 0;
  }
  return 1;
}

int main( int argc, char** argv ) {
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
  double* rspecg_initial = NULL;
  double rnorm[nfld] = {0., 0.};
  double expected_norm[nfld] = {0., 0.};
  double rspecg_norm[nfld] = {0., 0.};
  int owner[nfld] = {1, 1};

  if( trans.myproc == 1 ) {
    rspecg = calloc( nfld * trans.nspec2g, sizeof( double ) );
    rspecg_initial = calloc( nfld * trans.nspec2g, sizeof( double ) );
  }

  if( !has_zero_local_spectral && rspec == NULL ) {
    fprintf( stderr, "rank %d/%d: rspec allocation failed\n", trans.myproc, trans.nproc );
    return 1;
  }

  if( trans.myproc == 1 && rspecg == NULL ) {
    fprintf( stderr, "rank %d/%d: rspecg allocation failed\n", trans.myproc, trans.nproc );
    return 1;
  }

  if( trans.myproc == 1 && rspecg_initial == NULL ) {
    fprintf( stderr, "rank %d/%d: rspecg_initial allocation failed\n", trans.myproc, trans.nproc );
    return 1;
  }

  if( trans.myproc == 1 ) {
    setup_global_spectrum( rspecg, nfld, trans.nsmax );
    memcpy( rspecg_initial, rspecg, nfld * trans.nspec2g * sizeof( double ) );
    compute_global_spectral_norm_from_rspecg( rspecg, nfld, trans.nsmax, rspecg_norm );
    fprintf( stderr, "rank %d/%d: rspecg_norm = [", trans.myproc, trans.nproc );
    for( int field = 0; field < nfld; ++field )
      fprintf( stderr, " %g", rspecg_norm[field] );
    fprintf( stderr, " ]\n" );
  }
  compute_global_spectral_norm_from_spectral_values( nfld, trans.nsmax, expected_norm );

  if( trans.myproc == 1 ) {
    for( int field = 0; field < nfld; ++field ) {
      ASSERT( is_approx_eq( rspecg_norm[field], expected_norm[field] ) );
    }
  }

  fprintf( stderr,
           "rank %d/%d: nspec=%d nspec2=%d nspec2g=%d nspec2mx=%d nump=%d nfld=%d rspec=%p rspecg=%p\n",
           trans.myproc, trans.nproc, trans.nspec, trans.nspec2, trans.nspec2g, trans.nspec2mx,
           trans.nump, nfld, (void*) rspec, (void*) rspecg );

  if( trans.myproc == trans.nproc )
    ASSERT( has_zero_local_spectral );

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

  if( trans.myproc == 1 ) {
    for( int field = 0; field < nfld; ++field ) {
      ASSERT( is_approx_eq( rnorm[field], expected_norm[field] ) );
    }
  }

  fprintf( stderr, "rank %d/%d: trans_gathspec begin\n", trans.myproc, trans.nproc );
  struct GathSpec_t gathspec = new_gathspec( &trans );
    gathspec.rspec = rspec;
    gathspec.rspecg = rspecg;
    gathspec.nfld = nfld;
    gathspec.nto = owner;
  TRANS_CHECK( trans_gathspec( &gathspec ) );
  fprintf( stderr, "rank %d/%d: trans_gathspec end\n", trans.myproc, trans.nproc );

  if( trans.myproc == 1 ) {
    ASSERT( spectra_are_approx_eq( rspecg, rspecg_initial, nfld * trans.nspec2g ) );
    compute_global_spectral_norm_from_rspecg( rspecg, nfld, trans.nsmax, rspecg_norm );
    for( int field = 0; field < nfld; ++field ) {
      ASSERT( is_approx_eq( rspecg_norm[field], expected_norm[field] ) );
    }
  }

  free( rspec );
  free( rspecg );
  free( rspecg_initial );
  TRANS_CHECK( trans_delete( &trans ) );
  TRANS_CHECK( trans_finalize() );

  return 0;
}
