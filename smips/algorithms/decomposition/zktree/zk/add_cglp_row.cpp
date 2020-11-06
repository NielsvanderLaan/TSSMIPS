#include "zk.h"

void ZK::add_cglp_row(double *coef_x, double coef_theta, double *coef_y, double rhs)
{
  d_cglp.add_row(coef_x, coef_theta, coef_y, rhs);
}

void ZK::add_cglp_row(BendersCut &cut)
{
  vector<double> zeros(d_n2, 0.0);
  d_cglp.add_row(cut.d_beta.data(), 1 + cut.d_tau, zeros.data(), cut.d_alpha);
}