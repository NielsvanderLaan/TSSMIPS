#include "benders.h" 

BendersCut Benders::lpCut(double *x)
{ 
  double Tx[d_m2];
  computeTx(x, Tx);

  double dual[d_m2];
  double QLP = 0;
  fill(dual, dual + d_m2, 0.0);
    
  for (size_t s = 0; s != d_S; ++s)
  {
    double *ws = d_problem.d_omega[s].data(); // scenario (c-style array pointer)
    double prob = d_problem.d_probs[s]; 
    
    double rhs[d_m2];  // rhs vector (c-style array)
    for (size_t row = 0; row!= d_m2; ++row) // compute element-by-element
      rhs[row] = ws[row] - Tx[row];
    
    d_sub.update(rhs);
    Sub::Multipliers info = d_sub.solve();       
    
    double *lambda = info.lambda;
    QLP += prob * info.obj;
    
    for (size_t row = 0; row != d_m2; ++row)
      dual[row] += prob * lambda[row]; 

    delete[] lambda;
  }

  vector<double> beta(d_n1);
  for (size_t col = 0; col != d_n1; ++col)
  {
    for (size_t row = 0; row != d_m2; ++row)
      beta[col] += dual[row] * d_problem.d_Tmat[row][col];
  } 
    
  double alpha = inner_product(beta.begin(), beta.end(), x, QLP);

  return BendersCut{ alpha, beta, 0.0 };
}












