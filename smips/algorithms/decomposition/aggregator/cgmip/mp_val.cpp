#include "cgmip.h"

double CGMip::mp_val()
{
  if (mp_optimal())
    return d_mp.get(GRB_DoubleAttr_ObjVal);
  
  return 1e20;
}