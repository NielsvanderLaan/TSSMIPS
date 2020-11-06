#include "sub.h"

void Sub::update(vector<double> &x)
{
  vector<double> rhs = d_problem.d_omega[d_s];
  for (size_t row = 0; row != rhs.size(); ++row) // compute w - Tx element-by-element
    rhs[row] -= inner_product(d_problem.d_Tmat[row].begin(), d_problem.d_Tmat[row].end(), x.begin(), 0.0);

  d_model.set(GRB_DoubleAttr_RHS, d_constrs.data(), rhs.data(), rhs.size());

  d_model.update();
}