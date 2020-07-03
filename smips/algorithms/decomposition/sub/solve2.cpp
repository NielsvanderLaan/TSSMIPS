#include "sub.h"

Sub::GomInfo Sub::solve2()
{
  d_model.optimize();

  return GomInfo{d_model.get(GRB_DoubleAttr_Pi,  d_constrs.data(), d_m2),
                 d_model.get(GRB_IntAttr_VBasis, d_vars.data(),    d_n2),
                 d_model.get(GRB_IntAttr_CBasis, d_constrs.data(), d_m2)};
}