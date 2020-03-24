#include "cgmip.h"

bool CGMip::mp_optimal()
{
  return d_mp.get(GRB_IntAttr_Status) == 2;
}