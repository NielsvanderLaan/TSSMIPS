#include "benders.h" 

BendersCut Benders::sb_cut(double *x)
{
  double Tx[d_m2];
  computeTx(x, Tx);

  double alpha = 0.0;
  vector<double> beta(d_n1);
  
  double lw = 0;
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
       
    double pi[d_n1];  // pi = lambda T                      
    for (size_t var = 0; var != d_n1; ++var)
    {
      pi[var] = 0.0; 
      for (size_t row = 0; row != d_m2; ++row)
        pi[var] += lambda[row] * d_problem.d_Tmat[row][var];
      beta[var] += prob * pi[var];
    } 
    d_lr.update(ws, s, pi);
    double Lpiw = d_lr.solve();
    alpha += prob * Lpiw;



    

    for (size_t row = 0; row != d_m2; ++row)
      lw += prob * lambda[row] * ws[row]; 

    
    delete[] lambda;
  }

  return BendersCut{ alpha, beta, 0.0 };
}