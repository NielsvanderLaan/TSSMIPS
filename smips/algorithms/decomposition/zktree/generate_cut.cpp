#include "zktree.h"

BendersCut ZkTree::generate_cut(double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol)
{
  branch_cut(x, theta, rho, master, cuts, maxRounds, tol);

  for (size_t idx = 0; idx != d_nodes.size(); ++idx)
  {
    ZK* node = d_nodes[idx];
    if (not node->optimize()) // problem is infeasible?
    {
      cerr << "bac cut: subproblem is infeasible (CCR assumption not satisfied)\n";
      exit(1);
    }
    BendersCut cut = node->subgradient();
    add_row_to_cglp(cut.d_beta.data(), cut.d_tau, 1.0, cut.d_alpha, idx);
  }

  d_cglp.set(GRB_DoubleAttr_Obj, d_beta.data(), x, d_beta.size());
  d_tau.set(GRB_DoubleAttr_Obj, rho - d_L);

  solve_cglp();

  return d_candidate;
}