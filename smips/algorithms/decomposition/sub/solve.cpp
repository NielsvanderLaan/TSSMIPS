#include "sub.h"

Sub::Multipliers Sub::solve()
{
  d_model.optimize();
  return Multipliers {d_model.get(GRB_DoubleAttr_Pi, d_constrs, d_m2), d_model.get(GRB_DoubleAttr_ObjVal)};
}