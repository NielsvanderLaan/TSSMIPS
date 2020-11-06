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

  double *pi_ptr = d_model.get(GRB_DoubleAttr_Pi, d_constrs.data(), d_m2);
  vector<double> pi(pi_ptr, pi_ptr + d_m2);
  delete[] pi_ptr;

  return Multipliers { pi, d_model.get(GRB_DoubleAttr_ObjBound) };
}