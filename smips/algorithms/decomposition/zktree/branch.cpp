#include "zktree.h"

void ZkTree::branch(size_t node_idx, vector<double> &lb_nodes, double &UB, double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol)
{ 
  add_child(node_idx);   // appends child to d_nodes and updates cglp

  size_t var_idx = branch_var_idx(node_idx);     // uses the default strong branching
  
  size_t child_idx = d_nodes.size() - 1;
  
  update_bound(node_idx, var_idx, false);        // updates d_nodes[node_idx]: changes lower/upper bound of var (using floor and ceiling)                                    
  update_bound(child_idx, var_idx, true);        // also makes necessary changes to cglp
  
  update_global_bounds(node_idx, lb_nodes, UB, x, theta, rho, master, cuts, maxRounds, tol);
  update_global_bounds(child_idx, lb_nodes, UB, x, theta, rho, master, cuts, maxRounds, tol);
}    