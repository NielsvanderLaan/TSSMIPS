#include "cgmip.h"

CGMip::CGMip(const CGMip &other)
:
  d_mp(other.d_mp),
  d_beta(other.d_beta.size()),
  d_sub(other.d_sub),
  d_xVars(other.d_xVars.size()),
  d_yVars(other.d_yVars.size()),
  d_points(other.d_points)
{
  size_t n1 = d_beta.size();

  GRBVar *mp_vars = d_mp.getVars();
  d_alpha = mp_vars[0];
  copy_n(mp_vars + 1, n1, d_beta.begin());
  d_tau = mp_vars[n1 + 1]; 

  GRBVar *sub_vars = d_sub.getVars();
  copy_n(sub_vars, n1, d_xVars.begin());
  d_theta = sub_vars[n1];
  d_eta = sub_vars[n1 + 1];
  copy_n(sub_vars + n1 + 2, d_yVars.size(), d_yVars.begin());
  
  delete[] mp_vars;
  delete[] sub_vars; 
}