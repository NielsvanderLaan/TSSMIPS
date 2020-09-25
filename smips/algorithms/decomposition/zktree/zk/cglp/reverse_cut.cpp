#include "cglp.h"

void Cglp::reverse_cut(double UB)
{
  if (d_rcut_idx == -1)
  {
    d_rcut_idx = d_nMults;
    double *coeff_x = d_problem.d_c.data();

    for_each(coeff_x, coeff_x + d_n1, [](double &val){val *= -1;});
    vector<double> zeros(d_n2, 0.0);
    add_row(coeff_x, -1, zeros.data(), -UB);
    return;
  }

  size_t rhs_con = d_n1 + d_n2 + 1;
  double val = -UB + d_L;
  d_model.chgCoeff(d_constrs1[rhs_con], d_lambda1[d_rcut_idx], val);
  d_model.chgCoeff(d_constrs2[rhs_con], d_lambda2[d_rcut_idx], val);
  d_model.update();
}