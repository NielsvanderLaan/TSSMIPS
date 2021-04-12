#include "sub.h"

void Sub::update(vector<double> &x)
{
  vector<double> rhs = d_problem.d_omega[d_s];
  vector<vector<double>> &tech = d_problem.d_fix_tech ? d_problem.d_Tmat : d_problem.d_T_omega[d_s];
  for (size_t row = 0; row != rhs.size(); ++row) // compute w - Tx element-by-element
    rhs[row] -= inner_product(tech[row].begin(), tech[row].end(), x.begin(), 0.0);

  d_model.set(GRB_DoubleAttr_RHS, d_constrs.data(), rhs.data(), rhs.size());

  d_model.update();
}