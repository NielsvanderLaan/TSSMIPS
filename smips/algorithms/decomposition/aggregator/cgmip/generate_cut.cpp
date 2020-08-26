#include "cgmip.h"

BendersCut CGMip::generate_cut(double *x, double theta, bool init, double vwx, bool affine, double tol)
{
  d_tau.set(GRB_DoubleAttr_UB, affine ? 0 : GRB_INFINITY);      // force tau = 0 if only affine cuts are allowed
  
  set_mp_obj(x, theta);
  if (init)      // add initial point to mp to prevent unbounded rays
    add_mp_cut(Point{ vector<double>(x, x + d_xVars.size()), theta, vwx, 0.0, 0.0});

  BendersCut candidate{ 0, vector<double>(d_beta.size()), 0 };
  Point point{ vector<double>(d_xVars.size()), 0, 0, -GRB_INFINITY, GRB_INFINITY };

  while (true)
  {
    solve_mp();
    if (not mp_optimal())   // numerical issue occured
    {
      cout << "mp unbounded" << endl;
      break;
    }

    candidate = get_candidate();   // candidate cut

    set_sub_obj(candidate);        // attempt to find point which invalidates candidate cut
    point = solve_sub();
                                       // if cut is violated by point
    double diff = candidate.d_alpha - point.d_rhs_ub;
    if (diff > tol && check_mp_violation(max(diff - 1e-6, tol)))
    {
      add_mp_cut(point);           // add it to master
    }
    else   
      break;
  }


  candidate.d_alpha = point.d_rhs_lb;
  return candidate;       
}