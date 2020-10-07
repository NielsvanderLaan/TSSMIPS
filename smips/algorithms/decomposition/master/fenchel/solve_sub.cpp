#include "fenchel.h"

Fenchel::Point Fenchel::solve_sub()
{
  d_sub.optimize();

  double lb = d_sub.get(GRB_DoubleAttr_ObjBound);
  double ub;

  try
  {
    ub = d_sub.get(GRB_DoubleAttr_ObjVal);
  } catch (GRBException e)
  {
    cout << "error in Fenchel subproblem: did not find a feasible solution\n" << e.getErrorCode() << ' ' << e.getMessage() << endl;
    return Point {vector<double>(d_xvars.size()), 0, lb, GRB_INFINITY };
  }

  double *x = d_sub.get(GRB_DoubleAttr_X, d_xvars.data(), d_xvars.size());
  vector<double> xvals{ x, x + d_xvars.size() };
  delete[] x;
  return Point{ xvals, d_theta.get(GRB_DoubleAttr_X), lb, ub };
}