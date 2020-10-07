#include "tree.h"

bool Tree::solve(vector<Type> types, size_t node_idx, vector<double> &incumbent, double local_tol, double time_limit, bool rcuts, bool fenchel, size_t max_rounds)
{
  Benders *node = d_nodes[node_idx];
  double weight = 3.0;
  double upper_bound = (d_UB_global + weight * d_LB_nodes[node_idx]) / (1 + weight);    // * d_UB_global- global_tol
  Benders::Bounds bounds = node->hybrid_solve(types, true, max_rounds, upper_bound, local_tol, time_limit, rcuts, fenchel);

  d_LB_nodes[node_idx] = bounds.d_LB;

  d_LB_global = *min_element(d_LB_nodes.begin(), d_LB_nodes.end());
  if (bounds.d_UB < d_UB_global)
  {
    d_UB_global = bounds.d_UB;
    copy_n(node->d_incumbent, incumbent.size(), incumbent.begin());
    for (Benders *node: d_nodes)
      node->update(bounds.d_UB);
  }

  return bounds.branch;
}