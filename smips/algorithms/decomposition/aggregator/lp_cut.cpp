#include "aggregator.h"

BendersCut Aggregator::lp_cut(vector<double> &x)
{
  double dual[d_problem.d_m2];
  fill_n(dual, d_problem.d_m2, 0.0);
  double QLP = 0;

//#pragma omp parallel for reduction(+ : dual[:], QLP) num_threads(8)
  for (size_t s = 0; s < d_problem.d_S; ++s)
  {
    d_sub[s].update(x);
    Sub::Multipliers info = d_sub[s].solve();
    QLP += d_problem.d_probs[s] * info.obj;

    for (size_t row = 0; row != d_problem.d_m2; ++row)
      dual[row] += d_problem.d_probs[s] * info.lambda[row];
  }

  vector<double> beta(d_n1);
  for (size_t col = 0; col != d_n1; ++col)
  {
    for (size_t row = 0; row != d_problem.d_m2; ++row)
      beta[col] += dual[row] * d_problem.d_Tmat[row][col];
  }

  double alpha = inner_product(beta.begin(), beta.end(), x.begin(), QLP);

  return BendersCut{ alpha, beta, 0.0 };
}