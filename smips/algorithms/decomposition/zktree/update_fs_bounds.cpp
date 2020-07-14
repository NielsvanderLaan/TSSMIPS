#include "zktree.h"

void ZkTree::update_fs_bounds(size_t var, double val, bool lower)
{
  vector<vector<int>> &inds = lower ? d_lb_mult_inds : d_ub_mult_inds;
  for (size_t node_idx = 0; node_idx != d_nodes.size(); ++node_idx)
  {
    int mult_idx = inds[node_idx][var];
    if (mult_idx != -1)
    {
      d_cglp.chgCoeff(d_constrs[node_idx][d_beta.size() + 2], d_lambda[node_idx][mult_idx], val);
    } else
    {
      string name = "lambda_" + to_string(node_idx) + "_" + to_string(d_lambda[node_idx].size());
      GRBVar mult = d_cglp.addVar(lower ? 0 : -GRB_INFINITY, lower ? GRB_INFINITY : 0, 0, GRB_CONTINUOUS, name);
      d_lambda[node_idx].push_back(mult);
      inds[node_idx][var] = d_lambda[node_idx].size() - 1;
      d_cglp.chgCoeff(d_constrs[node_idx][d_beta.size() + 2], mult, val);
      d_cglp.chgCoeff(d_constrs[node_idx][var], mult, 1.0);
    }
  }
  d_cglp.update();

}