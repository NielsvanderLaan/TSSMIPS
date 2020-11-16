#include "zktree.h"

void ZkTree::optimize()
{
  d_cglp.optimize();

  if (d_cglp.get(GRB_IntAttr_Status) == 2)
  {
    d_cglp_val = d_cglp.get(GRB_DoubleAttr_ObjVal);

    double tau = d_tau.get(GRB_DoubleAttr_X);
    double *beta_ptr = d_cglp.get(GRB_DoubleAttr_X, d_beta.data(), d_beta.size());
    d_candidate = BendersCut {d_alpha.get(GRB_DoubleAttr_X) + tau * d_L + d_L, vector<double> (beta_ptr, beta_ptr + d_beta.size()), tau};
    delete[] beta_ptr;
  }
  else
  {
    d_cglp_val = GRB_INFINITY;
    cout << "error in ZkTree::optimize()\n";
    exit(d_cglp.get(GRB_IntAttr_Status));
  }
}