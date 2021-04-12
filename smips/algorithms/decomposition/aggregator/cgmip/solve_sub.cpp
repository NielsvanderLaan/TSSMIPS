#include "cgmip.h"

CGMip::Point CGMip::solve_sub(bool focus)
{
  d_sub.set(GRB_DoubleParam_MIPGap, focus ? 0.0 : 1e-4);

  d_sub.optimize();
  double rhs_lb = d_sub.get(GRB_DoubleAttr_ObjBound);
  double rhs_ub;

  try
  {
    rhs_ub = d_sub.get(GRB_DoubleAttr_ObjVal);
  } catch (GRBException e)
  {
    cout << "error in CGSP: did not find a feasible solution\n" << e.getErrorCode() << ' ' << e.getMessage() << endl;
    return Point {vector<double>(d_xVars.size()), 0, GRB_INFINITY, rhs_lb, GRB_INFINITY };
  }

  if (rhs_ub - rhs_lb > 1e-4 and not focus)
    return solve_sub(true);


  double *xVals = d_sub.get(GRB_DoubleAttr_X, d_xVars.data(), d_xVars.size());
  vector<double> x{ xVals, xVals + d_xVars.size() };
  delete[] xVals;


  return Point{ x, d_theta.get(GRB_DoubleAttr_X), d_eta.get(GRB_DoubleAttr_X), rhs_lb, rhs_ub };
}

