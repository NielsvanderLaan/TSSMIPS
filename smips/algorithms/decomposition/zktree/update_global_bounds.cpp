#include "zktree.h"

void ZkTree::update_global_bounds(size_t node_idx, vector<double> &lb_nodes, double &UB, double *x, double theta, Master &master, size_t maxRounds, bool gomory)
{
  if (node_idx >= lb_nodes.size())
    lb_nodes.push_back(-1);    

  if (solve(node_idx, x, theta, master, maxRounds, gomory)) // solve(node) returns true if LP-relaxation is feasible
  {
    ZK *node = d_nodes[node_idx];
    lb_nodes[node_idx] = node->d_objVal;     // update LB accordingly 
    if (is_feasible(node))                   // check if solution is integer
      UB = min(node->d_objVal, UB);          // update UB accordingly
  } else
    lb_nodes[node_idx] = GRB_INFINITY;               // LP-relaxation is infeasible
}