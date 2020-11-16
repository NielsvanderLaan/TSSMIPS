#include "zktree.h"

void ZkTree::update_global_bounds(size_t node_idx, vector<double> &lb_nodes, double &UB, double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol)
{
  if (node_idx >= lb_nodes.size())
    lb_nodes.push_back(-1);    

  if (solve(node_idx, x, theta, rho, master, cuts, maxRounds, tol)) // solve(node) returns true if LP-relaxation is feasible
  {
    ZK *node = d_nodes[node_idx];
    lb_nodes[node_idx] = node->d_objVal;     // update LB accordingly 
    if (is_feasible(node))                   // check if solution is integer
      UB = min(node->d_objVal, UB);          // update UB accordingly
  } else
    lb_nodes[node_idx] = GRB_INFINITY;       // LP-relaxation is infeasible
}