#include "zk.h"

double ZK::compute_a0(int row, vector<double> x, double theta)
{
  int e_i_ind[1] = {row};  double e_i_val[1] = {1.0};
  GRBsvec e_i {1, e_i_ind, e_i_val};  // unit vector
  int Brow_ind[d_nConstrs];  double Brow_val[d_nConstrs];
  GRBsvec Brow {d_nConstrs, Brow_ind, Brow_val};   // result vector
  GRBBSolve(d_model, &e_i, &Brow);   // extracting ith row of B^-1

  double rhs[d_nConstrs];

  for (size_t con = 0; con != d_nConstrs; ++con)
  {
    rhs[con] = d_omega[con] - d_tau[con] * theta;
    vector<double> &Trow = d_Tmat[con];
    for (size_t var = 0; var != Trow.size(); ++var)
      rhs[con] -= Trow[var] * x[var];
  }

  double ret = 0;

  for (size_t nz = 0; nz != Brow.len; ++nz)
    ret += Brow.val[nz] * rhs[Brow.ind[nz]];

  return ret;
}

