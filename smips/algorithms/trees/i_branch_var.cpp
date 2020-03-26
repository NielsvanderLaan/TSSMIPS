#include "tree.h"

int Tree::i_branch_var(double *x)
{
  size_t p1 = d_problem.d_p1;
  vector<double> scores(p1);
  bool integer = true;

  for (int var = 0; var != p1; ++var)
  {
    double val = x[var];
    if (not is_integer(val))
    {
      return var;
      integer = false;
      scores[var] = min(val - floor(val), ceil(val) - val);
    }
  }

  return integer ? -1 : distance(scores.begin(), max_element(scores.begin(), scores.end()));
}