#include "benderscallback.h"

void BendersCallback::callback()
{
  if (d_ncuts > 10)
    return;

  if (where == GRB_CB_MIPSOL)
    d_ben.reverse_cut(getDoubleInfo(GRB_CB_MIPSOL_OBJBST));

  if (where != GRB_CB_MIPNODE)
    return;

  if (getIntInfo(GRB_CB_MIPNODE_STATUS) != 2)
    return;

  double *xvals = getNodeRel(d_xvars.data(), d_xvars.size());
  vector<double> x(xvals, xvals + d_xvars.size());
  delete xvals;
  double theta = getNodeRel(d_theta);
  Master::Solution sol{x, theta, false};

  bool int_feas = all_of(x.begin(), x.begin() + d_problem.d_p1, [](double val){ return is_integer(val); });

  vector<double> vx = d_ben.d_agg.compute_vwx(x.data());
  double Qx = accumulate(vx.begin(), vx.end(),0.0) / d_problem.d_S;

  if (theta > Qx - 0.1)
    return;

  cout << "theta = " << theta << " Q(x) = " << Qx << '\n';

  BendersCut cut = d_ben.d_agg.strong_cut(sol, vx, false, 1e-4, int_feas);
  double cut_value = -inner_product(cut.d_beta.begin(), cut.d_beta.end(), x.begin(), -cut.d_alpha) / (1 + cut.d_tau);
  cout << "theta' = " << cut_value << '\n';


  d_ben.add_cut(cut, sol, 1e-4);

  GRBLinExpr betax;
  betax.addTerms(cut.d_beta.data(), d_xvars.data(), d_xvars.size());
  addLazy((1 + cut.d_tau) * d_theta + betax >= cut.d_alpha);


  ++d_ncuts;
}
