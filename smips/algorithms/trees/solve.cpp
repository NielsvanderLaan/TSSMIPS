#include "tree.h"

bool Tree::solve(size_t node_idx, vector<double> &incumbent, bool affine, double global_tol, double local_tol, int max_iter)
{
  Benders *node = d_nodes[node_idx];

  Benders::Bounds bounds = node->hybrid_solve(d_UB_global, affine, local_tol, max_iter);      // solve node
  
  if (bounds.d_infeasible)    // fathom if node is infeasible
  {                                          
    delete node;
    d_nodes.erase(d_nodes.begin() + node_idx);  
    d_LB_nodes.erase(d_LB_nodes.begin() + node_idx);
    return true;
  }
      // update LB, UB, and incumbent
  d_LB_nodes[node_idx] = bounds.d_LB;
  d_LB_global = *min_element(d_LB_nodes.begin(), d_LB_nodes.end());  
  if (bounds.d_UB < d_UB_global)
  {
    d_UB_global = bounds.d_UB;
    copy_n(node->d_incumbent, incumbent.size(), incumbent.begin());
  }

  return bounds.d_LB >= d_UB_global - global_tol;    // we should branch further: return true 
}