#include "tree.h"

bool Tree::solve(size_t node_idx, vector<double> &incumbent, bool affine, double global_tol, double local_tol)
{
  Benders *node = d_nodes[node_idx];
  Benders::Bounds bounds = node->hybrid_solve(d_UB_global, affine, d_problem.d_p1 > 0, local_tol);      // solve node

  d_LB_nodes[node_idx] = bounds.d_LB;

  d_LB_global = *min_element(d_LB_nodes.begin(), d_LB_nodes.end());
  if (bounds.d_UB < d_UB_global)
  {
    d_UB_global = bounds.d_UB;
    copy_n(node->d_incumbent, incumbent.size(), incumbent.begin());
  }

  if (bounds.d_LB > d_UB_global + 1e-8)    // fathom node (safe side here)
  {
    cout << "LB(node) = " << bounds.d_LB << '\n';
    delete node;
    d_nodes.erase(d_nodes.begin() + node_idx);
    d_LB_nodes.erase(d_LB_nodes.begin() + node_idx);
    return true;    // do not branch on this node
  }

  return bounds.d_LB > d_UB_global - global_tol;   // if true, then we do not branch further
}