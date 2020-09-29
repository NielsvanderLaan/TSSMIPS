#include "cglp.h"

void Cglp::reverse_cut(double UB)
{
  if (d_rcut_idx == -1)
  {
    d_rcut_idx = d_nMults;
    vector<double> zeros(d_n2, 0.0);
    add_row(d_problem.d_c.data(), 1, zeros.data(), UB, false);
    return;
  }

  size_t rhs_con = d_n1 + d_n2 + 1;
  double val = UB - d_L;
  d_model.chgCoeff(d_constrs1[rhs_con], d_lambda1[d_rcut_idx], val);
  d_model.chgCoeff(d_constrs2[rhs_con], d_lambda2[d_rcut_idx], val);
  d_model.update();

}