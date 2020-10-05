#include "cgmip.h"

void CGMip::update_bound(size_t var, double val, bool lower)
{
  if (lower)
    d_xVars[var].set(GRB_DoubleAttr_LB, val);
  else
    d_xVars[var].set(GRB_DoubleAttr_UB, val);

  GRBConstr *mp_cons = d_mp.getConstrs();
  for (size_t con = d_points.size() - 1; con != -1; --con)
  {
    Point &point = d_points[con];
    
    if ((lower && point.d_x[var] < val) || (not lower && point.d_x[var] > val))
    { 
      d_points.erase(d_points.begin() + con);
      d_mp.remove(mp_cons[con + 1]);
    }  
  }
  
  delete[] mp_cons;
  d_mp.update();
}