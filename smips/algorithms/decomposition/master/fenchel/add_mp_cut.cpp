#include "fenchel.h"

void Fenchel::add_mp_cut(const Fenchel::Point &point)
{
  d_points.push_back(point);

  GRBLinExpr lhs = point.d_theta * d_kappa;
  lhs.addTerms(point.d_xvals.data(), d_beta.data(), d_beta.size());
  d_mp.addConstr(lhs >= d_alpha);
  d_mp.update();
}