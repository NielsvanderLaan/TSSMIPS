#include "benders.h" 

BendersCut Benders::lbdaCut(double *x, double *alpha)
{
  double Tx[d_m2];
  computeTx(x, Tx);      // Tx is rba

      // cut coefficients: initialize to zero
  double gamma = 0.0;    
  double dual[d_m2];
  fill(dual, dual + d_m2, 0.0);

  for (size_t s = 0; s != d_S; ++s)
  {
    double *ws = d_problem.d_omega[s].data(); // scenario (c-style array pointer)
        // compute rhs, update subproblem
    
    double rhs[d_m2];          // rhs vector of subproblem (c-style array)  
    for (size_t row = 0; row != d_m2; ++row) // compute element-by-element
      rhs[row] = ws[row] - Tx[row];
    d_sub.update(rhs);           
    Sub::GomInfo info = d_sub.solve2();     // solve subproblem 
    
    double *lambda = info.lambda;           // extract lambda (for optimality cut)
    int *vBasis = info.vBasis;              // extract vBasis (to update gomory relaxation)   
    int *cBasis = info.cBasis;              // extract vBasis (to update gomory relaxation)     
    double gom_obj = compute_gomory(s, vBasis, cBasis, ws, alpha);
   
    double prob = d_problem.d_probs[s]; 
     
    gamma += prob * gom_obj;  // gom_obj = lambda^T (omega - alpha) + psi(omega - alpha), so we add lambda^T alpha in the following loop                              
  
    for (size_t row = 0; row != d_m2; ++row)
    {
      dual[row] += prob * lambda[row];  
      gamma += prob * lambda[row] * alpha[row];    
    }

    delete[] lambda; delete[] vBasis; delete[] cBasis;
  }

  vector<double> beta(d_n1);
  for (size_t col = 0; col != d_n1; ++col)
  {
    for (size_t row = 0; row != d_m2; ++row)
      beta[col] += dual[row] * d_problem.d_Tmat[row][col];
  } 
  
  return BendersCut{ gamma, beta, 0.0 };
}








