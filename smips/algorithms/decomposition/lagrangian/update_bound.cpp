#include "lagrangian.h"

void Lagrangian::update_bound(size_t var, double val, bool lower)
{
  if (lower)
    d_z_vars[var].set(GRB_DoubleAttr_LB, val);
  else
    d_z_vars[var].set(GRB_DoubleAttr_UB, val);

  d_model.update();
}