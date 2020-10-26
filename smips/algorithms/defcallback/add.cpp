#include "benderscallback.h"

void BendersCallback::add(BendersCut &cut, Master::Solution &sol, double tol)
{
  d_ben.add_cut(cut, sol, tol, true);

  GRBLinExpr betax;
  betax.addTerms(cut.d_beta.data(), d_xvars.data(), d_xvars.size());
  addLazy((1 + cut.d_tau) * d_theta + betax >= cut.d_alpha);
}