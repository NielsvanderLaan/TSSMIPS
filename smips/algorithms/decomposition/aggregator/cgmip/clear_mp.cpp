#include "cgmip.h"

void CGMip::clear_mp()
{
  GRBConstr *mp_cons = d_mp.getConstrs();
  for (size_t con = 0; con != d_points.size(); ++con)
    d_mp.remove(mp_cons[con + 1]);
  d_mp.update();
  delete[] mp_cons;

  d_points.clear();
}