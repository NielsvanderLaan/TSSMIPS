#include "cgmip.h"

void CGMip::clear_mp()
{
  GRBConstr *mp_cons = d_mp.getConstrs();
  for (size_t con = 0; con != d_points.size(); ++con)
    d_mp.remove(mp_cons[con]);
  d_mp.update();

  d_points.clear();
}