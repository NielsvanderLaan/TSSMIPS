#include "cgmip.h"

CGMip::Point CGMip::solve_sub()
{
  d_sub.optimize();
  
  int status = d_sub.get(GRB_IntAttr_Status);
  if (status != 2)
  {
    cout << "sub status: " << status << '\n';
    exit(1);
  }  
  
  double *xVals = d_sub.get(GRB_DoubleAttr_X, d_xVars.data(), d_xVars.size());
  vector<double> x{ xVals, xVals + d_xVars.size() };
  delete[] xVals;
  
  return Point{ x, d_theta.get(GRB_DoubleAttr_X), d_eta.get(GRB_DoubleAttr_X), d_sub.get(GRB_DoubleAttr_ObjBound), d_sub.get(GRB_DoubleAttr_ObjVal) };
}