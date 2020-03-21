#include "master.h"

BendersCut Master::transform_cut(double coef_theta, double *coef_x)
{
  double alpha = 1;
  double tau = coef_theta - 1;        
      // transforming cut (substituting slack variable expressions)
  for (size_t slack = 0; slack != d_nSlacks; ++slack)
  {
    double coeff = coef_x[d_n1 + slack];
    
    alpha += coeff * d_gamma[slack];
    tau += coeff * d_kappa[slack];
    
    for (size_t var = 0; var != d_n1; ++var)
      coef_x[var] += coeff * d_beta[slack][var];
  }
  
  vector<double> beta(coef_x, coef_x + d_n1);

  return BendersCut{ alpha, beta, tau };
}