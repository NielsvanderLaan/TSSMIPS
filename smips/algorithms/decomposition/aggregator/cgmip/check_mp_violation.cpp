#include "cgmip.h"

bool CGMip::check_mp_violation(double tol)
{
  if (d_mp.get(GRB_DoubleAttr_ConstrVio) > tol)
  {
    cout << "violation = " << d_mp.get(GRB_DoubleAttr_ConstrVio) << '\n';
    d_mp.set(GRB_IntParam_OutputFlag, 1);
    d_mp.optimize();
    exit(1);
  }

  return (d_mp.get(GRB_DoubleAttr_ConstrVio) <= tol);
}