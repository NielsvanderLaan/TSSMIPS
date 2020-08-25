#include "tree.h"

Tree::Split Tree::branch_var(size_t node_idx)
{
  Benders *node = d_nodes[node_idx];
  double *x = node->d_xvals;

  int var = i_branch_var(x);     // index of fractional-valued integer decision variable
  if (var != -1)
    return Split{ var, floor(x[var]), ceil(x[var]) };
  return Split{-1, -1, -1};   // no spatial branching


  var = c_branch_var(node, x);      // spatial branching variable idx
  if (var != -1)
  {
    if (var < d_problem.d_p1)
      return Split{var, ceil(x[var]) - 1, ceil(x[var])};

    return Split{var, x[var], x[var]};
  }

  /*
  var = c_branch_var_diam(node);
  if (var != -1)
  {
    double val = (node->d_lb[var] + node->d_ub[var]) / 2;

    if (var < d_problem.d_p1)
      return Split{var, ceil(val) - 1, ceil(val)};

    return Split{var, val, val};
  }
   */
  return Split{-1, -1, -1};
}