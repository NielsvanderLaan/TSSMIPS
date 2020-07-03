#include "zktree.h"

void ZkTree::update_fs_bounds(size_t var, double val, bool lower)
{
  vector<vector<int>> &inds = lower ? d_lb_mult_inds : d_ub_mult_inds;

  for (size_t node_idx = 0; node_idx = d_nodes.size(); ++node_idx)
  {
    for (size_t var = 0; var != d_beta.size(); ++var)
    {
      int mult_idx = inds[node_idx][var];
      if (mult_idx != -1)
      {
        d_cglp.chgCoeff(d_constrs[node_idx][d_beta.size() + 2], d_lambda[node_idx][mult_idx], val);
      } else
      {
        GRBVar mult = d_cglp.addVar(lower ? 0 : -GRB_INFINITY, lower ? GRB_INFINITY : 0, 0, GRB_CONTINUOUS);
        d_lambda[node_idx].push_back(mult);
        inds[node_idx][var] = d_lambda[node_idx].size() - 1;
        d_cglp.chgCoeff(d_constrs[node_idx][d_beta.size() + 2], mult, val);
        d_cglp.chgCoeff(d_constrs[node_idx][var], mult, val);
      }
    }
  }
}