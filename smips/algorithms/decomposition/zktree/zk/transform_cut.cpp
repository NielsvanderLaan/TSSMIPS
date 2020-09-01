#include "zk.h"

void ZK::transform_cut(double *coef_x, double *coef_y, double &coef_theta, double &coef_rhs, vector<double> &kappa, vector<vector<double>> &beta, vector<double> &gamma, size_t nSlacks)
{
  for (size_t slack = 0; slack != nSlacks; ++slack)
  {    
    double coef_slack = coef_x[d_n1 + slack];
    
    coef_theta += coef_slack * kappa[slack]; 
    
    for (size_t var = 0; var != d_n1; ++var)
      coef_x[var] += coef_slack * beta[slack][var]; 
   
    coef_rhs += coef_slack * gamma[slack];    // += because rhs is on the rhs      
    
    coef_x[d_n1 + slack] = 0;
  }

  transform_ycut(coef_x, coef_y, coef_theta, coef_rhs);
}