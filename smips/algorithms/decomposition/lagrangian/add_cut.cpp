#include "lagrangian.h"

void Lagrangian::add_cut(BendersCut &cut)
{
  GRBLinExpr lhs = (1 + cut.d_tau) * d_theta;
  lhs.addTerms(cut.d_beta.data(), d_z_vars.data(), d_z_vars.size());
  d_model->addConstr(lhs >= cut.d_alpha);
  d_model->update();
}