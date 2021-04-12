#include "sub.h"

vector<double> Sub::compute_slope(vector<double> &x)
{
  update(x);

  Multipliers info = solve();
  vector<double> &lambda = info.lambda;

  vector<double> pi(d_problem.d_n1, 0.0);  // pi = lambda T
  vector<vector<double>> &tech = d_problem.d_fix_tech ? d_problem.d_Tmat : d_problem.d_T_omega[d_s];

  for (size_t var = 0; var != pi.size(); ++var)
  {
    for (size_t row = 0; row != d_m2; ++row)
      pi[var] += lambda[row] * tech[row][var];
  }

  return pi;
}