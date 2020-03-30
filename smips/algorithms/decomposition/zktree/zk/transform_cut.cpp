#include "zk.h"

void ZK::transform_cut(double *coef_x, double *coef_y, double &coef_theta, double &coef_rhs, vector<double> &kappa, vector<vector<double>> &beta, vector<double> &gamma, size_t nSlacks)
{
  coef_rhs += coef_theta * d_L;
  for (size_t slack = 0; slack != nSlacks; ++slack)
  {    
    double coef_slack = coef_x[d_n1 + slack];
    
    coef_theta += coef_slack * kappa[slack]; 
    
    for (size_t var = 0; var != d_n1; ++var)
      coef_x[var] += coef_slack * beta[slack][var]; 
   
    coef_rhs += coef_slack * gamma[slack];    // += because rhs is on the rhs      
    
    coef_x[d_n1 + slack] = 0;
  }

  size_t slack = 0;
  for (size_t row = 0; row != d_nConstrs; ++row)
  {
    int sign = d_signs[row];
    if (sign == 0)      // equality constraint: no slack variable
      continue;
    
    double coeff = coef_y[d_n2 + slack];  // coefficient of slack variable  
      
    for (size_t var = 0; var != d_n1; ++var)
      coef_x[var] += sign * coeff * d_Tmat[row][var];
    for (size_t var = 0; var != d_n2; ++var)
      coef_y[var] += sign * coeff * d_Wmat[row][var];
      
    coef_theta += sign * coeff * d_tau[row];
    coef_rhs += sign * coeff * d_omega[row]; 
      
    ++slack; 
  }
}