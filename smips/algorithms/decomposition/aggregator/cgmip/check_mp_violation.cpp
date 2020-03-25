#include "cgmip.h"

bool CGMip::check_mp_violation(double tol)
{
  if (d_mp.get(GRB_DoubleAttr_ConstrVio) > tol)
  {
    cout << d_mp.get(GRB_DoubleAttr_ConstrVio) << '\n';
    exit(1);
  }

  return (d_mp.get(GRB_DoubleAttr_ConstrVio) <= tol);
}