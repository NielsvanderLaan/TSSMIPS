#include "fenchel.h"

void Fenchel::set_sub_obj(BendersCut &cut)
{
  d_theta.set(GRB_DoubleAttr_Obj, 1 + cut.d_tau);
  d_sub.set(GRB_DoubleAttr_Obj, d_xvars.data(), cut.d_beta.data(), d_xvars.size());

  d_sub.update();
}