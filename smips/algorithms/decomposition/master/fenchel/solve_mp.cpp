#include "fenchel.h"

bool Fenchel::solve_mp(double tol)
{
  d_mp.optimize();
  double violation = d_mp.get(GRB_DoubleAttr_ConstrVio) + d_mp.get(GRB_DoubleAttr_ConstrResidual);

  if (d_mp.get(GRB_IntAttr_Status) == 2 && (violation < tol))
    return true;

  chg_mp_tol(true);
  d_mp.optimize();
  chg_mp_tol(false);

  violation = d_mp.get(GRB_DoubleAttr_ConstrVio) + d_mp.get(GRB_DoubleAttr_ConstrResidual);

  return d_mp.get(GRB_IntAttr_Status) == 2 && (violation < tol);
}

