#include "cgmip.h"

void CGMip::add_row(BendersCut &cut)
{
  double kappa = 1 + cut.d_tau;
  GRBLinExpr lhs = kappa * d_theta;                        // kappa * theta + beta^T x
  lhs.addTerms(cut.d_beta.data(), d_xVars.data(), d_xVars.size());
  d_sub.addConstr(lhs, GRB_GREATER_EQUAL, cut.d_alpha);
  d_sub.update();

  if (cut.d_feas_cut) return;

  // check if feasibility cuts do not cut away points

  GRBConstr *constrs = d_mp.getConstrs();

  for (size_t con = 0; con != d_points.size(); ++con)
  {
    Point &point = d_points[con];
    double alpha_betax = -inner_product(point.d_x.begin(), point.d_x.end(), cut.d_beta.begin(), -cut.d_alpha);  // alpha - beta^T x
    point.d_theta = max(point.d_theta, alpha_betax / kappa);
    d_mp.chgCoeff(constrs[con], d_tau, point.d_theta);
  }

  delete[] constrs;
  d_mp.update();
}