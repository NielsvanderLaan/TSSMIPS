#include "cgmip.h"

bool CGMip::check_mp_violation(double tol)
{
  double violation = d_mp.get(GRB_DoubleAttr_ConstrVio) + d_mp.get(GRB_DoubleAttr_ConstrResidual);

  if (violation > tol)
  {
    if (not solve_mp(true))
    {
      cout << "mp unbounded\n";
      return false;
    }

    violation = d_mp.get(GRB_DoubleAttr_ConstrVio) + d_mp.get(GRB_DoubleAttr_ConstrResidual);
  }
  if (violation > tol)
    cout << "violation  = " << violation << endl;

  return (violation <= tol);
}