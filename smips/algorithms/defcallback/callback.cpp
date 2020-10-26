#include "benderscallback.h"

void BendersCallback::callback()
{
  if (where == GRB_CB_MIPSOL)
    d_ben.reverse_cut(getDoubleInfo(GRB_CB_MIPSOL_OBJBST));

  if (where != GRB_CB_MIPNODE)
    return;

  if (getIntInfo(GRB_CB_MIPNODE_STATUS) != 2)
    return;

  double tol = 1e-4;

  double *xvals = getNodeRel(d_xvars.data(), d_xvars.size());
  vector<double> x(xvals, xvals + d_xvars.size());
  delete xvals;
  double theta = getNodeRel(d_theta);
  Master::Solution sol{x, theta, false};

  bool int_feas = all_of(x.begin(), x.begin() + d_problem.d_p1, [](double val){ return is_integer(val); });

  vector<double> vx = d_ben.d_agg.compute_vwx(x.data());
  double Qx = accumulate(vx.begin(), vx.end(),0.0) / d_problem.d_S;

  if (theta > Qx - tol)
    return;

  cout << "Adding lazy cuts. theta = " << theta << " Q(x) = " << Qx << '\n';

  BendersCut cut;

  cut = d_ben.sb_cut(x.data());
  add(cut, sol, tol);

  cut = d_ben.d_agg.strong_cut(sol, vx, true, tol, int_feas);
  add(cut, sol, tol);

  cut = d_ben.d_agg.strong_cut(sol, vx, false, tol, int_feas);
  add(cut, sol, tol);
}
