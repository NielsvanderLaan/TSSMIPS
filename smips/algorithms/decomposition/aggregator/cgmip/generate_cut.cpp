#include "cgmip.h"

BendersCut CGMip::generate_cut(double *x, double theta, bool init, double vwx, bool affine, double tol, bool int_feas, double &gap)
{
  d_tau.set(GRB_DoubleAttr_UB, affine ? 0 : GRB_INFINITY);      // force tau = 0 if only affine cuts are allowed

  set_mp_obj(x, theta);
  if (init)      // add initial point to mp to prevent unbounded rays
    add_mp_cut(Point{ vector<double>(x, x + d_xVars.size()), theta, int_feas ? vwx : 1e6, 0.0, 0.0});

  BendersCut candidate{ 0, vector<double>(d_beta.size()), 0 };
  Point point{ vector<double>(d_xVars.size()), 0, GRB_INFINITY, -GRB_INFINITY, GRB_INFINITY };

  bool first_strike = false;
  while (true)
  {
    //d_mp.set(GRB_IntParam_OutputFlag, 1);
    if (not solve_mp(first_strike))
    {
      if (first_strike)
      {
        print("mp unbounded\n");
        break;
      }
      first_strike = true;
      continue;
    }
    candidate = get_candidate();   // candidate cut
    /*
    d_mp.write("mp.lp");
    GRBVar *vars = d_mp.getVars();
    for (size_t idx = 0; idx != d_mp.get(GRB_IntAttr_NumVars); ++ idx)
      cout << vars[idx].get(GRB_StringAttr_VarName) << ": " << vars[idx].get(GRB_DoubleAttr_X) << '\n';
    exit(123);
    */

    set_sub_obj(candidate);        // attempt to find point which invalidates candidate cut
    Point old_point = point;
    point = solve_sub();

    double diff = candidate.d_alpha - point.d_rhs_ub;
    if (diff < tol)     // optimal within tolerance
      break;
    if (distance(old_point, point) < 1e-8 || check_mp_violation(max(diff - 1e-6, tol)))
    {
      if (first_strike)
      {
        print("violation > tol\n");
        break;
      }
      first_strike = true;
      continue;
    }
    first_strike = false;
    add_mp_cut(point);
  }


  gap += candidate.d_alpha - point.d_rhs_lb;
  candidate.d_alpha = point.d_rhs_lb;
  return candidate;
}