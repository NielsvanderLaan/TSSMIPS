#include "cgmip.h"

void CGMip::set_sub_obj(BendersCut &cut)
{
  d_theta.set(GRB_DoubleAttr_Obj, cut.d_tau);
  d_sub.set(GRB_DoubleAttr_Obj, d_xVars.data(), cut.d_beta.data(), d_xVars.size());
}