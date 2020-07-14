#include "tree.h"

bool Tree::solve(size_t node_idx, vector<double> &incumbent, bool affine, double global_tol, double local_tol)
{
  Benders *node = d_nodes[node_idx];
  double weight = 3.0;
  double upper_bound = (d_UB_global + weight * d_LB_nodes[node_idx]) / (1 + weight);    // * d_UB_global- global_tol
  Benders::Bounds bounds = node->hybrid_solve(upper_bound, affine, true, local_tol);      // *


  d_LB_nodes[node_idx] = bounds.d_LB;

  cout << "LBs: ";
  for_each(d_LB_nodes.begin(), d_LB_nodes.end(), [](double val){ cout << val << ' '; });
  cout << '\n';

  d_LB_global = *min_element(d_LB_nodes.begin(), d_LB_nodes.end());
  if (bounds.d_UB < d_UB_global)
  {
    d_UB_global = bounds.d_UB;
    copy_n(node->d_incumbent, incumbent.size(), incumbent.begin());
  }


  return bounds.d_LB > d_LB_nodes[node_idx];  // *
  return bounds.d_LB > d_UB_global - global_tol;  // if true, then we do not branch further *
}