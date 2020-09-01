#include "zk.h"

bool ZK::gmi_cut(double *tab_row_x, double *tab_row_y, double a0, double *coef_x, double *coef_y, double &coef_theta, int nVarsMaster)
{
  double f0 = a0 - floor(a0);

  if (f0 < 1e-8  || 1 - f0 < 1e-8)
    return false;
      // theta
  coef_theta = max(tab_row_x[0] / f0, -tab_row_x[0] / (1 - f0));

      // x variables
  for (size_t var = 1; var != nVarsMaster; ++var)
  {
    size_t xvar = var - 1;
    
    double aj = tab_row_x[var];
    if (xvar < d_p1)    // check if integer variable
    {
      double fj = aj - floor(aj);
      coef_x[xvar] = min(fj / f0, (1 - fj) / (1 - f0));  
    } else 
    {
      coef_x[xvar]  = max(aj / f0, -aj / (1 - f0));
    }      
  } 
  
      // y variables
  for (size_t var = 0; var != d_nVars; ++var)
  {
    double aj = tab_row_y[var];
    
    if (var < d_p2)
    {
      double fj = aj - floor(aj);
      coef_y[var] = min(fj / f0, (1 - fj) / (1 - f0));  
    } else
    {
      coef_y[var] = max(aj / f0, -aj / (1 - f0));
    }
  }

  return true;
}