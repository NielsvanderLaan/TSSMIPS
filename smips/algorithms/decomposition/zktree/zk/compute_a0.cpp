#include "zk.h"

double ZK::compute_a0(int row, double yval, double theta, double rho)
{
  if (rho <= theta)
    return yval;

  int Brow_ind[d_nConstrs];  double Brow_val[d_nConstrs];
  GRBsvec Brow {d_nConstrs, Brow_ind, Brow_val};   // result vector
  B_inv(Brow, row);                             // extracting ith row of B^-1

  double Binv_tau = 0;
  for (int nz = 0; nz != Brow.len; ++nz)
    Binv_tau += Brow.val[nz] * d_tau[Brow.ind[nz]];

  return yval + Binv_tau * (rho - theta);

  /*
   *  Explanation: yval = {B^-1(omega - Tx - tau * rho)}_i (where i = row)
   *  We know yval, but we need
   *  a0 = {B^-1(omega - Tx - tau theta)}_i.
   *  To see this, see eq (2.2) in ZK and use the expression for d^k_i in Alg. 2, step 11, and use that
   *  A^-1 b = (x, theta)
   */
}

