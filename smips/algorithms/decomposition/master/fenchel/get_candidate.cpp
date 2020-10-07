#include "fenchel.h"

BendersCut Fenchel::get_candidate()
{
  double *b = d_mp.get(GRB_DoubleAttr_X, d_beta.data(), d_beta.size());
  vector<double> beta(b, b + d_beta.size());
  delete[] b;

  return BendersCut{ d_alpha.get(GRB_DoubleAttr_X), beta, d_kappa.get(GRB_DoubleAttr_X) - 1.0, true};
}