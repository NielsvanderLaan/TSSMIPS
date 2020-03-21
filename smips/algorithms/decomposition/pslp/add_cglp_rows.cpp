#include "pslp.h"

void Pslp::add_cglp_rows(double *coef_x, double coef_theta, double *coef_y, double rhs)
{
  for (size_t s = 0; s != d_S; ++s)
    d_zk[s].add_cglp_row(coef_x, coef_theta, coef_y, rhs);
}