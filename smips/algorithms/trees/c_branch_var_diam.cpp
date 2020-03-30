#include "tree.h"

int Tree::c_branch_var_diam(Benders *node)
{
  double best_score = -1;
  double best_idx = -1;

  for (int var = 0; var != d_problem.d_n1; ++var)
  {
    double score = (node->d_ub[var] - node->d_lb[var]) / (d_problem.d_u1[var] - d_problem.d_l1[var]);

    if (score > best_score)
    {
      best_score = score;
      best_idx = var;
    }
  }

  return best_idx;
}