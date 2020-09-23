#include "zktree.h"

void ZkTree::reverse_cut(double UB)     // reverse cut is of the form cx + theta <= UB
{
  if (d_rcut_inds[0] != -1)
  {
    for (size_t node_idx = 0; node_idx != d_nodes.size(); ++node_idx)
      d_cglp.chgCoeff(d_constrs[node_idx][d_beta.size() + 2], d_lambda[node_idx][d_rcut_inds[node_idx]], UB - d_L); // cx + theta' <= UB - d_L

  } else
  {
    for (size_t node_idx = 0; node_idx != d_nodes.size(); ++node_idx)
    {
      d_rcut_inds[node_idx] = d_lambda[node_idx].size();
      string name = "lambda_" + to_string(node_idx) + "_" + to_string(d_lambda[node_idx].size());
      GRBVar mult = d_cglp.addVar(-GRB_INFINITY, 0.0, 0, GRB_CONTINUOUS, name);   // <= constraint
      d_lambda[node_idx].push_back(mult);

      vector<GRBVar> lambda_vec(d_beta.size(), mult);
      d_cglp.chgCoeffs(d_constrs[node_idx].data(), lambda_vec.data(), d_problem.d_c.data(), d_beta.size());   // cx
      d_cglp.chgCoeff(d_constrs[node_idx][d_beta.size()], mult, 1.0);                                                   // + 1.0 * theta'
      d_cglp.chgCoeff(d_constrs[node_idx][d_beta.size() + 2], mult, UB -  d_L);                                         // <= UB - d_L
    }
  }

  d_cglp.update();
}

