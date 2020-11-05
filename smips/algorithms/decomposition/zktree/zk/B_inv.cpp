#include "zk.h"

void ZK::B_inv(GRBsvec &result, int row)
{
  int e_i_ind[1] = {row};
  double e_i_val[1] = {1.0};
  GRBsvec e_i {1, e_i_ind, e_i_val};  // unit vector

  GRBBSolve(d_model, &e_i, &result);   // extracting ith row of B^-1
}