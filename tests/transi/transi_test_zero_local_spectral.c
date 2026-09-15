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

#ifndef LAM_VERSION
static int setup_sh_global_spectrum( double* rspecg, int nfld, int nsmax ) {
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
  return coeff;
}

static void compute_sh_global_spectral_norm_from_spectral_values( int nfld, int nsmax, double* rnorm ) {
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

static void compute_sh_global_spectral_norm_from_rspecg( const double* rspecg, int nfld, int nsmax, double* rnorm ) {
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

#else

static int lam_nmax( int tx, int ty, int m ) {
  return (int) ( (double) ty / tx * sqrt( (double) ( tx * tx - m * m ) ) + 1.e-10 );
}

static int setup_lam_global_spectrum( double* rspecg, int nfld, int tx, int ty ) {
  int coeff = 0;
  for( int m = 0; m <= tx; ++m ) {
    for( int n = 0; n <= lam_nmax( tx, ty, m ); ++n ) {
      for( int part = 0; part < 4; ++part ) {
        for( int field = 0; field < nfld; ++field ) {
          rspecg[spectral_index( coeff, field, nfld )] = spectral_value( field, m, n, part );
        }
        ++coeff;
      }
    }
  }
  return coeff;
}

static void compute_lam_global_spectral_norm_from_spectral_values( int nfld, int tx, int ty, double* rnorm ) {
  for( int field = 0; field < nfld; ++field ) {
    rnorm[field] = 0.;
    for( int m = 0; m <= tx; ++m ) {
      for( int n = 0; n <= lam_nmax( tx, ty, m ); ++n ) {
        for( int part = 0; part < 4; ++part ) {
          const double value = spectral_value( field, m, n, part );
          rnorm[field] += value * value;
        }
      }
    }
    rnorm[field] = sqrt( rnorm[field] );
  }
}

static void compute_lam_global_spectral_norm_from_rspecg(
    const double* rspecg, int nfld, int tx, int ty, double* rnorm ) {
  for( int field = 0; field < nfld; ++field ) {
    rnorm[field] = 0.;
  }

  int coeff = 0;
  for( int m = 0; m <= tx; ++m ) {
    for( int n = 0; n <= lam_nmax( tx, ty, m ); ++n ) {
      for( int part = 0; part < 4; ++part ) {
        (void) n;
        (void) part;
        for( int field = 0; field < nfld; ++field ) {
          const double value = rspecg[spectral_index( coeff, field, nfld )];
          rnorm[field] += value * value;
        }
        ++coeff;
      }
    }
  }

  for( int field = 0; field < nfld; ++field ) {
    rnorm[field] = sqrt( rnorm[field] );
  }
}
#endif

static int setup_global_spectrum( double* rspecg, int nfld, struct Trans_t* trans ) {
#ifndef LAM_VERSION
  return setup_sh_global_spectrum( rspecg, nfld, trans->nsmax );
#else
  return setup_lam_global_spectrum( rspecg, nfld, trans->nmsmax, trans->nsmax );
#endif
}

static void compute_global_spectral_norm_from_spectral_values( int nfld, struct Trans_t* trans, double* rnorm ) {
#ifndef LAM_VERSION
  compute_sh_global_spectral_norm_from_spectral_values( nfld, trans->nsmax, rnorm );
#else
  compute_lam_global_spectral_norm_from_spectral_values( nfld, trans->nmsmax, trans->nsmax, rnorm );
#endif
}

static void compute_global_spectral_norm_from_rspecg( const double* rspecg, int nfld, struct Trans_t* trans, double* rnorm ) {
#ifndef LAM_VERSION
  compute_sh_global_spectral_norm_from_rspecg( rspecg, nfld, trans->nsmax, rnorm );
#else
  compute_lam_global_spectral_norm_from_rspecg( rspecg, nfld, trans->nmsmax, trans->nsmax, rnorm );
#endif
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
#ifndef LAM_VERSION
  TRANS_CHECK( trans_set_trunc( &trans, 2 ) );
#else
  TRANS_CHECK( trans_set_trunc_lam( &trans, 2, 3 ) );
#endif
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
    ASSERT( setup_global_spectrum( rspecg, nfld, &trans ) == trans.nspec2g );
    memcpy( rspecg_initial, rspecg, nfld * trans.nspec2g * sizeof( double ) );
    compute_global_spectral_norm_from_rspecg( rspecg, nfld, &trans, rspecg_norm );
    fprintf( stderr, "rank %d/%d: rspecg_norm = [", trans.myproc, trans.nproc );
    for( int field = 0; field < nfld; ++field )
      fprintf( stderr, " %g", rspecg_norm[field] );
    fprintf( stderr, " ]\n" );
  }
  compute_global_spectral_norm_from_spectral_values( nfld, &trans, expected_norm );

  if( trans.myproc == 1 ) {
    for( int field = 0; field < nfld; ++field ) {
      ASSERT( is_approx_eq( rspecg_norm[field], expected_norm[field] ) );
    }
  }

  fprintf( stderr,
           "rank %d/%d: nspec=%d nspec2=%d nspec2g=%d nspec2mx=%d nump=%d nfld=%d rspec=%p rspecg=%p\n",
           trans.myproc, trans.nproc, trans.nspec, trans.nspec2, trans.nspec2g, trans.nspec2mx,
           trans.nump, nfld, (void*) rspec, (void*) rspecg );

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
      ASSERT( isfinite( rnorm[field] ) );
      ASSERT( rnorm[field] > 0. );
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
    compute_global_spectral_norm_from_rspecg( rspecg, nfld, &trans, rspecg_norm );
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
