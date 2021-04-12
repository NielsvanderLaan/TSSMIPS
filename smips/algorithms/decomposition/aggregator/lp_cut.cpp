#include "aggregator.h"

BendersCut Aggregator::lp_cut(vector<double> &x)
{

  double beta[d_n1];
  fill_n(beta, d_n1, 0.0);
  double QLP = 0;

#pragma omp parallel for reduction(+ : beta[:], QLP)
  for (size_t s = 0; s < d_problem.d_S; ++s)
  {
    d_sub[s].update(x);
    Sub::Multipliers info = d_sub[s].solve();
    QLP += d_problem.d_probs[s] * info.obj;

    double prob = d_problem.d_probs[s];
    vector<vector<double>> &tech = d_problem.d_fix_tech ? d_problem.d_Tmat : d_problem.d_T_omega[s];
    for (size_t col = 0; col != d_n1; ++col)
    {
      for (size_t row = 0; row != d_problem.d_m2; ++row)
        beta[col] += prob * info.lambda[row] * tech[row][col];
    }


  }

  double alpha = inner_product(beta, beta + d_n1, x.begin(), QLP);
  return BendersCut{ alpha, vector<double>{beta, beta + d_n1}, 0.0 };
}