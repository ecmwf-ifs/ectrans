/*
 * (C) Copyright 2026- ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ectrans/transi.h"
#include "transi_test.h"

// -------------------------------------------------------------------------------------------------

void test_trans_inquire_lam() {
  // -----------------------------------------------------------------------------------------------
  // Test local transforms
  // -----------------------------------------------------------------------------------------------

  struct Trans_t trans;

  const int nx = 20;
  const int ny = 18;
  const double dx = 2500.0;
  const double dy = 2500.0;
  const int tx = (nx - 1) / 2;
  const int ty = (ny - 1) / 2;

  TRANS_CHECK(trans_new(&trans));
  TRANS_CHECK(trans_set_resol_lam(&trans, nx, ny, dx, dy));
  TRANS_CHECK(trans_set_trunc_lam(&trans, ty, tx));
  TRANS_CHECK(trans_setup(&trans));

  // Check all keys are accepted
  TRANS_CHECK(trans_inquire(&trans, "numpp,ngptotl,nmyms,npossp,nptrms,nallms,ndim0g,nvalue,mvalue,"
    "nfrstlat,nlstlat,nptrlat,nptrfrstlat,nptrlstlat,nsta,nonl,nultpp,nptrls,npms,ndgl,nsmax,"
    "myproc,nproc,llam,nspec,nspec2,nspec2g,nspec2mx,nump,ngptot,ngptotg,ngptotmx,n_regions_ns,"
    "n_regions_ew,my_region_ns,my_region_ew,nfrstloff,nptrfloff,nprtrns,nlei3,nspolegl,nmsmax"));

  TRANS_CHECK(trans_delete(&trans));
}

int main(int argc, char **argv) {
  trans_use_mpi( test_use_mpi() );

  setbuf(stdout,NULL); // unbuffered stdout

  printf("-----------------------------\n");
  test_trans_inquire_lam();

  TRANS_CHECK(trans_finalize());

  return 0;
}
