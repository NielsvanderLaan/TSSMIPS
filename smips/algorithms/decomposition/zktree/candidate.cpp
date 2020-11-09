#include "zktree.h"

BendersCut ZkTree::candidate()
{
  if (d_cglp.get(GRB_IntAttr_Status) != 2)
  {
    cerr << "error in cglp of zktree, zee ZkTree::candidate()\n";
    exit(1);
  }

  double alpha = d_alpha.get(GRB_DoubleAttr_X);
  double tau = d_tau.get(GRB_DoubleAttr_X);

  double *beta_ptr = d_cglp.get(GRB_DoubleAttr_X, d_beta.data(), d_beta.size());
  vector<double> beta(beta_ptr, beta_ptr + d_beta.size());
  delete[] beta_ptr;

  return BendersCut {alpha + tau * d_L + d_L, beta, tau};
}