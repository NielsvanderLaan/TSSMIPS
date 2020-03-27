#include "cgmip.h"

void CGMip::remove_mp_cuts()
{
  d_mp.optimize();
  
  if (not mp_optimal()) return;

  GRBConstr *cons = d_mp.getConstrs();
  for (size_t row = d_points.size() - 1; row != - 1; --row)   
  {
    if (cons[row].get(GRB_IntAttr_CBasis) == 0)
    {
      d_mp.remove(cons[row]);
      d_points.erase(d_points.begin() + row);  
    }   
  }
  
  d_mp.update();
  delete[] cons;
}