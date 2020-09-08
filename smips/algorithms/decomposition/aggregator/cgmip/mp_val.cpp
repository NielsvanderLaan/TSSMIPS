#include "cgmip.h"

double CGMip::mp_val()
{
  if (mp_optimal())
    return d_mp.get(GRB_DoubleAttr_ObjVal) +d_alpha.get(GRB_DoubleAttr_X) - d_sub.get(GRB_DoubleAttr_ObjBound) ;

  return GRB_INFINITY;
}