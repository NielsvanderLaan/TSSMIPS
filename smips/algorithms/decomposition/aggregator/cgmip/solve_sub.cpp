#include "cgmip.h"

CGMip::Point CGMip::solve_sub()
{
  d_sub.optimize();
  
  int status = d_sub.get(GRB_IntAttr_Status);
  if (status != 2)
  {
    cout << "sub status: " << status << '\n';
    d_sub.write("sub.lp");
    exit(1);
  }

  double rhs_lb = d_sub.get(GRB_DoubleAttr_ObjBound);
  double rhs_ub = d_sub.get(GRB_DoubleAttr_ObjVal);
  if (rhs_lb < rhs_ub - 1e-4)
  {
    d_sub.set(GRB_DoubleParam_MIPGap, 0.0);
    d_sub.optimize();
    d_sub.set(GRB_DoubleParam_MIPGap, 1e-4);
    rhs_lb = d_sub.get(GRB_DoubleAttr_ObjBound);
    rhs_ub = d_sub.get(GRB_DoubleAttr_ObjVal);
  }
  
  double *xVals = d_sub.get(GRB_DoubleAttr_X, d_xVars.data(), d_xVars.size());
  vector<double> x{ xVals, xVals + d_xVars.size() };
  delete[] xVals;



  return Point{ x, d_theta.get(GRB_DoubleAttr_X), d_eta.get(GRB_DoubleAttr_X), rhs_lb, rhs_ub };
}