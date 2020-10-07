#include "fenchel.h"

BendersCut Fenchel::fenchel_cut(vector<double> &x, double theta, double tol)
{
  cout << "Fenchel::fenchel_cut()" << endl;
  set_mp_obj(x, theta);

  BendersCut candidate{ 0, vector<double>(d_beta.size()), 0 };
  Point point{ vector<double>(d_xvars.size()), 0, -GRB_INFINITY, GRB_INFINITY };

  while (true)
  {
    if (not solve_mp())     // mp status is not optimal
      break;

    candidate = get_candidate();
    set_sub_obj(candidate);
    point = solve_sub();

    // TODO: guard against looping (identifying same point and candidate over and over)
    // store previous candidate and compute difference + check mp tolerance.
    // clear mp if numerical issues occur

    if (candidate.d_alpha - point.d_ub < tol)
      break;
    else
      add_mp_cut(point);
  }

  candidate.d_alpha = point.d_lb;
  cout << "done" << endl;
  return candidate;
}