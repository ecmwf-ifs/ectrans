// (C) Copyright 2026- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.

#include "hicblas.h"

extern "C" {
void mult_butm_sp(
  char transpose, const int *order, const int *levels, const int *betalen_max,
  const int *lev_offset, const int *lev_ij, const int *lev_ik, const int *lev_ibetalen,
  const int *lev_node_offset, const int *lev_node_ifcol, const int *lev_node_ilcol,
  const int *lev_node_ifrow, const int *lev_node_ilrow, const int *lev_node_icols,
  const int *lev_node_irows, const int *lev_node_irank, const int *lev_node_ioffbeta,
  const int *lev_node_iclist_offset, const int *lev_node_iclist, const int *lev_node_pnonim_offset,
  const float *lev_node_pnonim, const int *lev_node_b_offset, const float *lev_node_b,
  const float *A, float *C) {
}

void mult_butm_dp(
  char transpose, const int *order, const int *levels, const int *betalen_max,
  const int *lev_offset, const int *lev_ij, const int *lev_ik, const int *lev_ibetalen,
  const int *lev_node_offset, const int *lev_node_ifcol, const int *lev_node_ilcol,
  const int *lev_node_ifrow, const int *lev_node_ilrow, const int *lev_node_icols,
  const int *lev_node_irows, const int *lev_node_irank, const int *lev_node_ioffbeta,
  const int *lev_node_iclist_offset, const int *lev_node_iclist, const int *lev_node_pnonim_offset,
  const double *lev_node_pnonim, const int *lev_node_b_offset, const double *lev_node_b,
  const double *A, double *C) {
}
}
