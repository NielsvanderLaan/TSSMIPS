#include "benders.h" 

BendersCut Benders::lbdaCut(vector<double> &x, vector<double> &alpha)
{
      // cut coefficients: initialize to zero
  double gamma = 0.0;    
  double dual[d_m2];
  fill(dual, dual + d_m2, 0.0);

  for (size_t s = 0; s != d_S; ++s)
  {
    Sub &sub = d_agg.d_sub[s];
    sub.update(x);
    Sub::GomInfo info = sub.solve2();

    double gom_obj = compute_gomory(s, info.vbasis, info.cbasis, d_problem.d_omega[s], alpha);

    gamma += d_problem.d_probs[s] * inner_product(info.lambda.begin(), info.lambda.end(), alpha.begin(), gom_obj);

    for (size_t row = 0; row != d_m2; ++row)
      dual[row] += d_problem.d_probs[s] * info.lambda[row];
  }

  vector<double> beta(d_n1);
  for (size_t col = 0; col != d_n1; ++col)
  {
    for (size_t row = 0; row != d_m2; ++row)
      beta[col] += dual[row] * d_problem.d_Tmat[row][col];
  } 
  
  return BendersCut{ gamma, beta, 0.0 };
}








