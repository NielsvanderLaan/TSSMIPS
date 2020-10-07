#include "fenchel.h"

void Fenchel::add_row(BendersCut &cut)
{
  if (cut.d_feas_cut) return;

  double kappa = 1 + cut.d_tau;
  GRBLinExpr lhs = kappa * d_theta;                        // kappa * theta + beta^T x
  lhs.addTerms(cut.d_beta.data(), d_xvars.data(), d_xvars.size());
  d_sub.addConstr(lhs, GRB_GREATER_EQUAL, cut.d_alpha);

  GRBConstr *mp_cons = d_mp.getConstrs();
  for (size_t con = 0; con != d_points.size(); ++con)
  {
    Point &point = d_points[con];
    double alpha_betax = -inner_product(point.d_xvals.begin(), point.d_xvals.end(), cut.d_beta.begin(), -cut.d_alpha);  // alpha - beta^T x
    point.d_theta = max(point.d_theta, alpha_betax / kappa);
    d_mp.chgCoeff(mp_cons[con], d_kappa, point.d_theta);
  }
  delete[] mp_cons;

  d_sub.update();
  d_mp.update();
}