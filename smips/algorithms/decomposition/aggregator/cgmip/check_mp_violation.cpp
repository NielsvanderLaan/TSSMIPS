#include "cgmip.h"

bool CGMip::check_mp_violation(double tol)
{
  return (d_mp.get(GRB_DoubleAttr_ConstrVio) <= tol);
}