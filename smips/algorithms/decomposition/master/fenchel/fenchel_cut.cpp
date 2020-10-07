#include "fenchel.h"

BendersCut Fenchel::fenchel_cut(vector<double> &x, double theta, double tol, bool reset)
{
  if (reset) clear_mp();
  set_mp_obj(x, theta);

  BendersCut candidate{ 0, vector<double>(d_beta.size()), 0 };
  Point point{ vector<double>(d_xvars.size()), 0, -GRB_INFINITY, GRB_INFINITY };

  while (true)
  {
    if (not solve_mp(tol))     // numerical issues
    {
      if (reset)
        break;
      return fenchel_cut(x, theta, tol, true);
    }

    candidate = get_candidate();
    set_sub_obj(candidate);
    point = solve_sub(tol);

    if (candidate.d_alpha - point.d_ub < tol)
      break;
    else
      add_mp_cut(point);
  }

  candidate.d_alpha = point.d_lb;
  return candidate;
}