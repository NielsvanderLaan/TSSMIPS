#include "sub.h"

void Sub::update(double *rhs)
{
  d_model.set(GRB_DoubleAttr_RHS, d_constrs.data(), rhs, d_m2);
}