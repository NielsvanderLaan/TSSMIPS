#include "cgmip.h"

static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}

double CGMip::mp_max_coeff()
{
  double *beta = d_mp.get(GRB_DoubleAttr_X, d_beta.data(), d_beta.size());
  double beta_max = abs(*max_element(beta , beta + d_beta.size(), abs_compare));
  delete[] beta;

  double alpha = abs(d_alpha.get(GRB_DoubleAttr_X));
  double tau = abs(d_tau.get(GRB_DoubleAttr_X));

  return max(max(alpha, tau), beta_max);
}