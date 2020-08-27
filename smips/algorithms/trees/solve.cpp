#include "tree.h"

bool Tree::solve(bool lp_cuts, bool sb_cuts, bool zk_cuts, bool strong_cuts,
                 bool affine, size_t node_idx, vector<double> &incumbent, double local_tol)
{
  Benders *node = d_nodes[node_idx];
  double weight = 3.0;
  double upper_bound = (d_UB_global + weight * d_LB_nodes[node_idx]) / (1 + weight);    // * d_UB_global- global_tol
  Benders::Bounds bounds = node->hybrid_solve(lp_cuts, sb_cuts, zk_cuts, strong_cuts, upper_bound, affine, local_tol);

  d_LB_nodes[node_idx] = bounds.d_LB;


  d_LB_global = *min_element(d_LB_nodes.begin(), d_LB_nodes.end());
  if (bounds.d_UB < d_UB_global)
  {
    d_UB_global = bounds.d_UB;
    copy_n(node->d_incumbent, incumbent.size(), incumbent.begin());
  }


  return bounds.branch;
}