#include "cgmip.h"

BendersCut CGMip::generate_cut(double *x, double theta, bool init, double vwx, bool affine, double tol, bool int_feas, double &gap)
{
  d_tau.set(GRB_DoubleAttr_UB, affine ? 0 : GRB_INFINITY);      // force tau = 0 if only affine cuts are allowed

  set_mp_obj(x, theta);
  if (init)      // add initial point to mp to prevent unbounded rays
    add_mp_cut(Point{ vector<double>(x, x + d_xVars.size()), theta, int_feas ? vwx : 1e7, 0.0, 0.0});

  BendersCut candidate{ 0, vector<double>(d_beta.size()), 0 };
  Point point{ vector<double>(d_xVars.size()), 0, GRB_INFINITY, -GRB_INFINITY, GRB_INFINITY };

  bool first_strike = false;
  while (true)
  {
    if (not solve_mp(first_strike))
    {
      if (first_strike)
      {
        cout << "mp unbounded\n";
        break;
      }
      first_strike = true;
      continue;
    }
    candidate = get_candidate();   // candidate cut

    set_sub_obj(candidate);        // attempt to find point which invalidates candidate cut
    Point old_point = point;
    point = solve_sub();

    if (candidate.d_alpha - point.d_rhs_ub < tol)     // optimal within tolerance
      break;
    if (distance(old_point, point) < 1e-8)
    {
      if (first_strike)
      {
        cout << "violation > tol\n";
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