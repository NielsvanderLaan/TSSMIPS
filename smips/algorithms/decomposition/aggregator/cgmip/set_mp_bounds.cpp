#include "cgmip.h"

void CGMip::set_mp_bounds(double M)
{
  vector<double> lb (d_beta.size(), -M);
  vector<double> ub (d_beta.size(), M);

  d_mp.set(GRB_DoubleAttr_LB, d_beta.data(), lb.data(), lb.size());
  d_mp.set(GRB_DoubleAttr_UB, d_beta.data(), ub.data(), ub.size());

  d_alpha.set(GRB_DoubleAttr_LB, -M);
  d_alpha.set(GRB_DoubleAttr_UB, M);

  d_tau.set(GRB_DoubleAttr_LB, 0);
  d_tau.set(GRB_DoubleAttr_UB, M);
}