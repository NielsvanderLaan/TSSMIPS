#include "zktree.h"

size_t ZkTree::branch_var_idx(size_t node_idx, bool strong)
{
  ZK *node = d_nodes[node_idx];
  vector<double> &yvals = node->d_yvals;
  
  if (not strong)
  {
    for (size_t var = 0; var != node->d_p2; ++var)
    {
      if (not is_integer(yvals[var]))
        return var;
    }
    return -1;  
  }
  
  
  
  double eps = 1e-6;
  double best_score = -1;
  double best_idx = -1;
  
  for (size_t var = 0; var != node->d_p2; ++var)
  {
    if (not is_integer(yvals[var]))
    {
      double delta_plus = node->probe(var, ceil(yvals[var]), true);      // increase in LP-objective value by imposing y[var] >= ceil(yvals[var])
      double delta_minus = node->probe(var, floor(yvals[var]), false);   // increase in LP-objective value by imposing y[var] <= floor(yvals[var]) 
      if (delta_plus + delta_minus > 1e10)    // one of the children is infeasible
        return var;                           // branch on this variable
      
      double score = max(delta_minus, eps) * max(delta_plus, eps);      // compute score (using rule proposed in PhD thesis Tobias Achterberg)
      if (score > best_score)    // compare to best score
      {
        best_score = score;      // if better, update best score  
        best_idx = var;          // and corresponding index  
      }
    }
  }
  
  return best_idx;
}