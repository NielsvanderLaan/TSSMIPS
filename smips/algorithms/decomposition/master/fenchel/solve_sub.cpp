#include "fenchel.h"

Fenchel::Point Fenchel::solve_sub(double tol, bool focus)
{
  d_sub.set(GRB_DoubleParam_MIPGap, focus ? 0 : 1e-4);
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

  if (ub - lb > tol and not focus)
    return solve_sub(tol, true);

  double *x = d_sub.get(GRB_DoubleAttr_X, d_xvars.data(), d_xvars.size());
  vector<double> xvals{ x, x + d_xvars.size() };
  delete[] x;
  return Point{ xvals, d_theta.get(GRB_DoubleAttr_X), lb, ub };
}