#include "zk.h"

void ZK::add_cglp_row(double *coef_x, double coef_theta, double *coef_y, double rhs)
{
  d_cglp.add_row(coef_x, coef_theta, coef_y, rhs);
}