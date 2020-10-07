#include "master.h"

void Master::add_cut(BendersCut &cut)
{
  d_fenchel.add_row(cut);
  double kappa = 1 + cut.d_tau;
  ++d_nSlacks;

  // adding the cut to d_cmodel
  GRBaddvar(d_cmodel, 0, NULL, NULL, 0, 0, GRB_INFINITY, GRB_CONTINUOUS, NULL);  // slack

  size_t numVars = d_n1 + 2; // theta, x-vars, and slack

  int cind[numVars];
  iota(cind, cind + d_n1 + 1, 0);
  cind[d_n1 + 1] = d_n1 + d_nSlacks;

  double cval[numVars];
  cval[0] = kappa;
  copy_n(cut.d_beta.begin(), d_n1, cval + 1);

  cval[numVars - 1] = -1;     // >= constraint, so slack features with -1
  GRBaddconstr(d_cmodel, numVars, cind, cval, GRB_EQUAL, cut.d_alpha - kappa * d_L, NULL);

  // update slack identities
  d_kappa.push_back(kappa);
  d_beta.push_back(cut.d_beta);
  d_gamma.push_back(cut.d_alpha);

  GRBupdatemodel(d_cmodel);
}