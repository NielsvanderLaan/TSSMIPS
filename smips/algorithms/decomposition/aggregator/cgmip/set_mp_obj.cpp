#include "cgmip.h"

void CGMip::set_mp_obj(double *x, double &theta)
{
  d_mp.set(GRB_DoubleAttr_Obj, d_beta.data(), x, d_beta.size());
  d_tau.set(GRB_DoubleAttr_Obj, theta);
  d_mp.update();
}