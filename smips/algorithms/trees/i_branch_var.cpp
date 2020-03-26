#include "tree.h"

int Tree::i_branch_var(double *x)
{
  for (int var = 0; var != d_problem.d_p1; ++var)
  {
    if (not is_integer(x[var]))
      return var;
  }

  return -1;
}