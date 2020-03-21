#include "cgmip.h"

void CGMip::add_mp_cut(Point const &point)
{
  // add the cut -alpha + x^T beta + theta tau >= -eta 
  GRBLinExpr lhs = -d_alpha + point.d_theta * d_tau;
  lhs.addTerms(point.d_x.data(), d_beta.data(), point.d_x.size());
  d_mp.addConstr(lhs, GRB_GREATER_EQUAL, -point.d_eta);
  d_mp.update();
  d_points.push_back(point);
}