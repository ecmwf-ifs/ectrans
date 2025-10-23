#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ectrans/transi.h"

#define TRANS_CHECK( CALL ) do {\
  int errcode = CALL;\
  if( errcode != TRANS_SUCCESS) {\
    printf("ERROR: %s failed @%s:%d:\n%s\n",#CALL,__FILE__,__LINE__,trans_error_msg(errcode));\
    exit(1);\
  }\
} while(0)

const int truncation = 79;

int main() {
  // Get handle to instance of ecTrans
  struct Trans_t trans;
  TRANS_CHECK(trans_new(&trans));

  // Calculate parameters of full Gaussian grid
  int nlat = 2 * (truncation + 1);
  int nlon = 2 * nlat;

  // Set resolution of ecTrans
  int* nloen = malloc(sizeof(int) * nlat);
  for (int i = 0; i < nlat; ++i) {
    nloen[i] = nlon;
  }
  TRANS_CHECK(trans_set_resol(&trans, nlat, nloen));
  TRANS_CHECK(trans_set_trunc(&trans, truncation));

  // Set up ecTrans
  TRANS_CHECK(trans_setup(&trans));

  printf("nspec2 = %d\n", trans.nspec2);
  printf("ngptot = %d\n", trans.ngptot);

  // Allocate our work arrays
  double* spectral_field = calloc(trans.nspec2, sizeof(double));
  double* spectral_field_2 = calloc(trans.nspec2, sizeof(double));
  double* grid_point_field = malloc(trans.ngptot * sizeof(double));

  // Initialise our spectral array
  int m = 3, n = 8;
  TRANS_CHECK(trans_inquire(&trans, "nasm0")); // Populate the "nasm0" member of trans
  spectral_field[trans.nasm0[m] + 2 * (n - m) + 1] = 1.0;

  // Perform an inverse transform
  struct InvTrans_t invtrans = new_invtrans(&trans);
  invtrans.nscalar = 1;
  invtrans.rspscalar = spectral_field;
  invtrans.rgp = grid_point_field;
  trans_invtrans(&invtrans);

  // Perform a direct transform back
  struct DirTrans_t dirtrans = new_dirtrans(&trans);
  dirtrans.nscalar = 1;
  dirtrans.rgp = grid_point_field;
  dirtrans.rspscalar = spectral_field_2;
  trans_dirtrans(&dirtrans);

  // Compute error between first and second spectral field arrays
  double norm2 = 0.0;
  for (int i = 0; i < trans.nspec2; ++i) {
    norm2 += pow(spectral_field_2[i] - spectral_field[i], 2);
  }
  norm2 = sqrt(norm2);
  printf("Error = %e\n", norm2);

  // Deallocate arrays
  free(grid_point_field);
  free(spectral_field_2);
  free(spectral_field);
  free(nloen);

  // Free ecTrans
  TRANS_CHECK(trans_delete(&trans));

  return 0;
}
