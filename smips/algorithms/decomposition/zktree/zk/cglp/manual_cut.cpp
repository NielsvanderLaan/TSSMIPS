#include "cglp.h"

Cut Cglp::manual_cut(size_t var_idx)
{
  size_t nl1 = d_lambda1.size();
  size_t nl2 = d_lambda2.size();
  
  double *lambda1 = d_model.get(GRB_DoubleAttr_X, d_lambda1.data(), nl1);
  double *lb1 = d_model.get(GRB_DoubleAttr_LB, d_lambda1.data(), nl1);
  double *ub1 = d_model.get(GRB_DoubleAttr_UB, d_lambda1.data(), nl1);

  transform(lambda1, lambda1 + nl1, lb1, lambda1, [](double a, double b){ return max(a, b); });
  transform(lambda1, lambda1 + nl1, ub1, lambda1, [](double a, double b){ return min(a, b); });

  double *lambda2 = d_model.get(GRB_DoubleAttr_X, d_lambda2.data(), nl2);
  double *lb2 = d_model.get(GRB_DoubleAttr_LB, d_lambda2.data(), nl2);
  double *ub2 = d_model.get(GRB_DoubleAttr_UB, d_lambda2.data(), nl2);

  transform(lambda2, lambda2 + nl2, lb2, lambda2, [](double a, double b){ return max(a, b); });
  transform(lambda2, lambda2 + nl2, ub2, lambda2, [](double a, double b){ return min(a, b); });

  size_t dim = d_constrs1.size() - 1;
  double coeff[dim];
  vector<double> lambda1A(dim, 0.0), lambda2A(dim, 0.0);

  double rhs;
  for (size_t con = 0; con != dim; ++con)
  {
    for (size_t var = 0; var != nl1; ++var)
      lambda1A[con] += d_model.getCoeff(d_constrs1[con], d_lambda1[var]) * lambda1[var];
    for (size_t var = 0; var != nl2; ++var)
      lambda2A[con] += d_model.getCoeff(d_constrs2[con], d_lambda2[var]) * lambda2[var];
    coeff[con] = max(lambda1A[con], lambda2A[con]);
  }

  double left = 0;
  double right = 0;
  for (size_t var = 0; var != nl1; ++var)
    left += d_model.getCoeff(d_constrs1.back(), d_lambda1[var]) * lambda1[var];
  for (size_t var = 0; var != nl2; ++var)
    right += d_model.getCoeff(d_constrs2.back(), d_lambda2[var]) * lambda2[var];

  rhs = min(left, right);

  vector<double> Trow(coeff, coeff + d_n1);
  double r = coeff[d_n1];
  vector<double> Wrow(coeff + d_n1 + 1, coeff + d_n1 + 1 + d_n2);
  rhs += r*d_L;

    // strengthening
  double lambda1_0 = lambda1[0];
  double lambda2_0 = lambda2[0];

  double diff = lambda1_0 - lambda2_0;

  if (diff > -1e-8)
    return Cut{ Trow, r, Wrow, rhs};

  for (size_t var = 0; var != d_p1; ++var)
  {
    double m = (lambda2A[var] - lambda1A[var]) / diff;
    Trow[var] = min(lambda1A[var] + lambda1_0 * floor(m), lambda2A[var] + lambda2_0 * ceil(m));
  }

  for (size_t var = 0; var != d_p2; ++var)
  {
    if (var == var_idx)
      continue;

    size_t con = d_n1 + 1 + var;   // first d_n1 + 1 constraints correspond to x and theta, we need the y variables
    double m = (lambda2A[con] - lambda1A[con]) / diff;
    Wrow[var] = min(lambda1A[con] + lambda1_0 * floor(m), lambda2A[con] + lambda2_0 * ceil(m));
  }

  delete[] lambda1;
  delete[] lb1;
  delete[] ub1;
  delete[] lambda2;
  delete[] lb2;
  delete[] ub2;

  return Cut {Trow, r, Wrow, rhs};
}