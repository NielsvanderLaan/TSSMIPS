#include "cgmip.h"

BendersCut CGMip::generate_cut(double *x, double theta, double vwx, bool affine, double tol, bool int_feas, double &gap, bool reset)
{
  if (reset) clear_mp();
  d_tau.set(GRB_DoubleAttr_UB, affine ? 0 : GRB_INFINITY);      // force tau = 0 if only affine cuts are allowed

  set_mp_obj(x, theta);
  d_S.set(GRB_DoubleAttr_LB, int_feas ? -vwx : -1e6);

  BendersCut candidate{ 0, vector<double>(d_beta.size()), 0 };
  Point point{ vector<double>(d_xVars.size()), 0, GRB_INFINITY, -GRB_INFINITY, GRB_INFINITY };

  bool first_strike = false;
  while (true)
  {
    if (not solve_mp(first_strike))
    {
      print("mp unbounded: resolving with more focus\n");
      if (not first_strike)
      {
        first_strike = true;
        continue;
      }
      if (not reset)
      {
        print("mp unbounded: resetting\n");
        return generate_cut(x, theta, vwx, affine, tol, int_feas, gap, true);
      }
      print("mp unbounded (after reset)");
      break;
    }

    candidate = get_candidate();   // candidate cut

    set_sub_obj(candidate);        // attempt to find point which invalidates candidate cut
    Point old_point = point;
    point = solve_sub();

    double diff = candidate.d_alpha - point.d_rhs_ub;
    //cout << "diff = " << diff << '\n';
    if (diff < tol)     // optimal within tolerance
      break;
    if (distance(old_point, point) < 1e-8 || check_mp_violation(max(diff - 1e-6, tol)))
    {
      if (not first_strike)
      {
        print("violation > tol: resolving with more focus\n");
        first_strike = true;
        continue;
      }
      if (not reset)
      {
        print("violation > tol: resetting\n");
        return generate_cut(x, theta, vwx, affine, tol, int_feas, gap, true);
      }
      print("violation > tol (after reset)\n");
      break;
    }

    first_strike = false;
    add_mp_cut(point);
  }

  gap += candidate.d_alpha - point.d_rhs_lb;
  candidate.d_alpha = point.d_rhs_lb;

  return candidate;
}