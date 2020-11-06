#include "lagrangian.h"

void Lagrangian::update(vector<double> &pi)
{
  d_model.set(GRB_DoubleAttr_Obj, d_z_vars.data(), pi.data(), pi.size());
  d_model.update();
}