#include "fenchel.h"

void Fenchel::update_bounds(int var, double val, bool lower)
{
  d_xvars[var].set(lower ? GRB_DoubleAttr_LB : GRB_DoubleAttr_UB, val);

  GRBConstr *mp_cons = d_mp.getConstrs();
  for (size_t idx = 0; idx != d_points.size(); ++idx)
  {
    Point &point = d_points[idx];
    if ((lower && point.d_xvals[var] < val) || (not lower && point.d_xvals[var] > val))
    {
      d_points.erase(d_points.begin() + idx);
      d_mp.remove(mp_cons[idx]);
    }
  }

  delete[] mp_cons;

  d_mp.update();
  d_sub.update();
}