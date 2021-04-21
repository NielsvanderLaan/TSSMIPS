#include "cgmip.h"

BendersCut CGMip::generate_cut(double *x, double theta, bool init, double vwx, bool affine, double tol, bool int_feas, double &gap,
                               double &niter, double &sptime, double &mptime,
                               bool reset)
{
  niter = 0; sptime = 0; mptime = 0;
  if (reset) clear_mp();
  d_tau.set(GRB_DoubleAttr_UB, affine ? 0 : GRB_INFINITY);      // force tau = 0 if only affine cuts are allowed

  set_mp_obj(x, theta);
  if (init)      // add initial point to mp to prevent unbounded rays
    add_mp_cut(Point{ vector<double>(x, x + d_xVars.size()), theta, int_feas ? vwx : 1e6, 0.0, 0.0});

  BendersCut candidate{ 0, vector<double>(d_beta.size()), 0 };
  Point point{ vector<double>(d_xVars.size()), 0, GRB_INFINITY, -GRB_INFINITY, GRB_INFINITY };

  bool first_strike = false;
  while (true)
  {
    ++niter;
    auto mp_t1 = chrono::high_resolution_clock::now();
    bool mp_succes = solve_mp(first_strike, affine);
    auto mp_t2 = chrono::high_resolution_clock::now();
    mptime += chrono::duration_cast<chrono::milliseconds>(mp_t2 - mp_t1).count() / 1000.0;

    if (not mp_succes)
    {
      if (not first_strike)
      {
        //print("mp unbounded: resolving mp with focus\n");
        first_strike = true;
        continue;
      }
      if (not reset)
      {
        print("mp unbounded: resetting\n");
        return generate_cut(x, theta, true, vwx, affine, tol, int_feas, gap, niter, sptime, mptime, true);
      }

      print("mp unbounded (after reset)\n");
      break;
    }

    candidate = get_candidate();   // candidate cut

    set_sub_obj(candidate);        // attempt to find point which invalidates candidate cut
    Point old_point = point;

    auto sp_t1 = chrono::high_resolution_clock::now();
    point = solve_sub();
    auto sp_t2 = chrono::high_resolution_clock::now();
    sptime += chrono::duration_cast<chrono::milliseconds>(sp_t2 - sp_t1).count() / 1000.0;

    double diff = candidate.d_alpha - point.d_rhs_ub;
    if (diff < tol)     // optimal within tolerance
      break;

    if (distance(old_point, point) < 1e-8|| check_mp_violation(max(diff - 1e-6, tol)))
    {
      if (not first_strike)
      {
        //print("violation > tol: resolving mp with focus\n");
        first_strike = true;
        continue;
      }
      if (not reset)
      {
        print("violation > tol: resetting\n");
        return generate_cut(x, theta, true, vwx, affine, tol, int_feas, gap, niter, sptime, mptime, true);
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