#include "cgmip.h"

bool CGMip::check_mp_violation(double tol)
{
  double violation = d_mp.get(GRB_DoubleAttr_ConstrVio) + d_mp.get(GRB_DoubleAttr_ConstrResidual);

  if (violation > tol)
    cout << "violation = " << violation << '\n';

  return (violation <= tol);
}