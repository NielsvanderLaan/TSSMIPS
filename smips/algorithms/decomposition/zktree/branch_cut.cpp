#include "zktree.h"

void ZkTree::branch_cut(double *x, double theta, Master &master, size_t maxRounds, bool gomory)
{
  double tol = 1e-4;  
  double UB = GRB_INFINITY;
  
  vector<double> lb_nodes(d_nodes.size());      // objective value of LP-relaxation of node  
  
  for (size_t node_idx = 0; node_idx != d_nodes.size(); ++node_idx)
  {
    d_nodes[node_idx]->update(x, theta);
    update_global_bounds(node_idx, lb_nodes, UB, x, theta, master, maxRounds, gomory);
  }
  
  auto node_iterator = min_element(lb_nodes.begin(), lb_nodes.end());
  size_t node_idx = distance(lb_nodes.begin(), node_iterator);
  double LB = *node_iterator;    
  
  while (LB < UB - tol)
  { 
    if (lb_nodes[node_idx] >= UB)    // sub-optimal node, skip this one (but NEVER fathom)
      continue; 
    
    branch(node_idx, lb_nodes, UB, x, theta, master, maxRounds, gomory);    
    
    node_iterator = min_element(lb_nodes.begin(), lb_nodes.end());
    node_idx = distance(lb_nodes.begin(), node_iterator);
    LB = *node_iterator;
  }
  
  cout << "LB = " << LB << ". UB = " << UB << ".\n";
  cout << "number of nodes: " << d_nodes.size() << '\n';
}













