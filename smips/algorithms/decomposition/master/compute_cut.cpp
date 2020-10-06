#include "master.h"

void Master::compute_cut(vector<double> &tab_row, double a0, double &coef_theta, double *coef_x)
{
  double a_theta = tab_row[0];

  double f0 = a0 - floor(a0);

  coef_theta = max(a_theta / f0, -a_theta / (1 - f0));
  
  for (size_t var = 1; var != tab_row.size(); ++var)
  {
    size_t xvar = var - 1;
    double aj = tab_row[var];
    if (xvar < d_p1)    // integer variable
    {
      double fj = aj - floor(aj);
      coef_x[xvar] = min(fj / f0, (1 - fj) / (1 - f0));  
    } else               // continuous variable
      coef_x[xvar] = max(aj / f0, -aj / (1 - f0));   
  }
}