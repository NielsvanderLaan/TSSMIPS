#include "tree.h"

Tree::Split Tree::branch_var(size_t node_idx)
{
  Benders *node = d_nodes[node_idx];
  double *x = node->d_xvals;

  int var = i_branch_var(x);     // index of fractional-valued integer decision variable    
  if (var != -1)
    return Split{ var, floor(x[var]), ceil(x[var]) };

  var = c_branch_var(node, x);      // spatial branching variable idx
  if (var == -1)
    var = c_branch_var_diam(node);

  if (var == -1)
    return Split { -1, -1, -1 };

  double val = x[var];
  if (var < d_problem.d_p1)      // spatial branching on integer variable
    return Split{var, floor(val), floor(val) + 1};

  return Split{var, val, val};
}