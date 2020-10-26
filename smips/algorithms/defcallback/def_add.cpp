#include "DEF.h"

void DEF::add(BendersCut &cut)
{
  GRBLinExpr lhs = (1 + cut.d_tau) * d_theta;
  lhs.addTerms(cut.d_beta.data(), d_xvars.data(), d_xvars.size());

  d_model.addConstr(lhs >= cut.d_alpha);
  d_model.update();
}