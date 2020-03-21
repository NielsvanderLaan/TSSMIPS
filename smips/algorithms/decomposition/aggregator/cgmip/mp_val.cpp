#include "cgmip.h"

double CGMip::mp_val()
{
  return d_mp.get(GRB_DoubleAttr_ObjVal);
}