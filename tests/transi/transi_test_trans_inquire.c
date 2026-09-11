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

void test_trans_inquire() {
  // -----------------------------------------------------------------------------------------------
  // Test global transforms
  // -----------------------------------------------------------------------------------------------

  struct Trans_t trans;

  const unsigned int truncation = 79;
  const unsigned int nlat = 2 * (truncation + 1);
  const unsigned int nlon = 2 * nlat;
  int* nloen  = malloc(sizeof(int) * nlat);
  for (int i = 0; i < nlat; ++i) {
    nloen[i] = nlon;
  }

  TRANS_CHECK(trans_new(&trans));
  TRANS_CHECK(trans_set_resol(&trans, nlat, nloen));
  TRANS_CHECK(trans_set_trunc(&trans, truncation));
  TRANS_CHECK(trans_setup(&trans));

  // Check all keys are accepted
  TRANS_CHECK(trans_inquire(&trans, "numpp,ngptotl,nmyms,nasm0,npossp,nptrms,nallms,ndim0g,nvalue,"
    "nfrstlat,nlstlat,nptrlat,nptrfrstlat,nptrlstlat,nsta,nonl,nultpp,nptrls,nnmeng,rmu,rgw,rpnm,"
    "npms,rlapin,ndglu"));

  TRANS_CHECK(trans_delete(&trans));
  free(nloen);
}

int main(int argc, char **argv) {
  trans_use_mpi( test_use_mpi() );

  setbuf(stdout,NULL); // unbuffered stdout

  printf("-----------------------------\n");
  test_trans_inquire();

  TRANS_CHECK(trans_finalize());

  return 0;
}
