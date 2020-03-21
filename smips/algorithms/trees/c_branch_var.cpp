#include "tree.h"

int Tree::c_branch_var(Benders *node, double *x)
{ 
  double best_score = -1;
  double best_idx = -1;


  for (int var = 0; var != d_problem.d_n1; ++var)
  {
    if (x[var] == node->d_lb[var] || x[var] == node->d_ub[var])
      continue;
    
    double lb = d_problem.d_l1[var];
    double ub = d_problem.d_u1[var];
    double score = min(x[var] - lb, ub - x[var]) / (ub - lb);
    
    if (score > best_score)
    {
      best_score = score;
      best_idx = var;
    }
  }
  
  return best_idx;
}