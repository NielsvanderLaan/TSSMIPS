#include "master.h"

void Master::strengthen_cut(BendersCut &cut)
{
  d_interceptor.set(GRB_DoubleAttr_Obj, d_xvars.data(), cut.d_beta.data(), d_xvars.size());
  d_theta.set(GRB_DoubleAttr_Obj, 1 + cut.d_tau);

  d_interceptor.optimize();
  double alpha_prime = d_interceptor.get(GRB_DoubleAttr_ObjBound);


  cut.d_alpha = max(cut.d_alpha, alpha_prime);

}