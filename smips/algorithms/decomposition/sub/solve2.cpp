#include "sub.h"

Sub::GomInfo Sub::solve2()
{
  d_model.optimize();

  double *pi_ptr = d_model.get(GRB_DoubleAttr_Pi, d_constrs.data(), d_m2);
  vector<double> pi(pi_ptr, pi_ptr + d_m2);
  delete[] pi_ptr;

  int *v_ptr = d_model.get(GRB_IntAttr_VBasis, d_vars.data(),  d_n2);
  vector<int> vbasis(v_ptr, v_ptr + d_n2);
  delete[] v_ptr;

  int *c_ptr = d_model.get(GRB_IntAttr_CBasis, d_constrs.data(),  d_m2);
  vector<int> cbasis(c_ptr, c_ptr + d_m2);
  delete[] c_ptr;

  return GomInfo{pi, vbasis, cbasis};
}