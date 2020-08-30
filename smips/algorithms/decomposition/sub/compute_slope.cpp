#include "sub.h"

vector<double> Sub::compute_slope(size_t s, double *x)
{
  vector<double> rhs = d_problem.d_omega[s];
  for (size_t row = 0; row != rhs.size(); ++row) // compute w - Tx element-by-element
    rhs[row] -= inner_product(d_problem.d_Tmat[row].begin(), d_problem.d_Tmat[row].end(), x, 0.0);

  update(rhs.data(), s);
  Multipliers info = solve();
  double *lambda = info.lambda;

  vector<double> pi(d_problem.d_n1, 0.0);  // pi = lambda T
  for (size_t var = 0; var != pi.size(); ++var)
  {
    for (size_t row = 0; row != d_m2; ++row)
      pi[var] += lambda[row] * d_problem.d_Tmat[row][var];
  }

  delete[] lambda;
  return pi;
}