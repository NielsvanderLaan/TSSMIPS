#include "sub.h"

Sub::Multipliers Sub::solve()
{
  d_model.optimize();
  int status = d_model.get(GRB_IntAttr_Status);
  if (status == 3)
  {
    cout << "subproblem is infeasible: no complete recourse. Implement feasibility cuts\n";
    exit(1);
  }

  return Multipliers {d_model.get(GRB_DoubleAttr_Pi, d_constrs, d_m2), d_model.get(GRB_DoubleAttr_ObjVal)};
}