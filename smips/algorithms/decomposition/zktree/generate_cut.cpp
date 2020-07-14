#include "ZkTree.h"

BendersCut ZkTree::generate_cut(double *x, double theta, Master &master, size_t maxRounds, bool gomory)
{
  branch_cut(x, theta, master, maxRounds, gomory);

  for (size_t idx = 0; idx != d_nodes.size(); ++idx)
  {
    ZK* node = d_nodes[idx];
    node->optimize(); // what if problem is infeasible?
    BendersCut cut = node->subgradient();
    //cout << -inner_product(cut.d_beta.begin(), cut.d_beta.end(), x, -cut.d_alpha) << '\n';
    add_row_to_cglp(cut.d_beta.data(), cut.d_tau, 1.0, cut.d_alpha, idx);
  }

  d_cglp.set(GRB_DoubleAttr_Obj, d_beta.data(), x, d_beta.size());
  d_tau.set(GRB_DoubleAttr_Obj, theta);

  d_cglp.optimize();
  //cout << "CGLP = " << d_cglp.get(GRB_DoubleAttr_ObjVal) << '\n';

  double *beta_ptr = d_cglp.get(GRB_DoubleAttr_X, d_beta.data(), d_beta.size());
  vector<double> beta(beta_ptr, beta_ptr + d_beta.size());
  delete[] beta_ptr;


  double tau = d_tau.get(GRB_DoubleAttr_X);
  double alpha = d_alpha.get(GRB_DoubleAttr_X) + tau * d_L + d_L;

  return BendersCut { alpha, beta, tau };
}