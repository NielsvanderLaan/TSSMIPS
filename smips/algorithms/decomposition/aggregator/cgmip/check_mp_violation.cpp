#include "cgmip.h"

bool CGMip::check_mp_violation(double tol)
{
  double violation = d_mp.get(GRB_DoubleAttr_ConstrVio) + d_mp.get(GRB_DoubleAttr_ConstrResidual);


  if (violation > tol)
  {
    d_mp.reset();
    d_mp.set(GRB_IntParam_ScaleFlag, 0);
    d_mp.set(GRB_IntParam_NumericFocus, 3);
    solve_mp();
    d_mp.set(GRB_IntParam_ScaleFlag, -1);
    d_mp.set(GRB_IntParam_NumericFocus, 0);
    violation = d_mp.get(GRB_DoubleAttr_ConstrVio) + d_mp.get(GRB_DoubleAttr_ConstrResidual);
    if (violation > tol)
      cout << "violation  = " << violation << endl;
  }


  return (violation <= tol);
}