#include "sub.h"

Sub::Sub(const Sub &other):
d_model(other.d_model),
d_problem(other.d_problem),
d_s(other.d_s),
d_m2(other.d_m2),
d_n2(other.d_n2)
{
  GRBVar *var_ptr = d_model.getVars();
  for (size_t var = 0; var != d_n2; ++var)
    d_vars.push_back(var_ptr[var]);
  delete[] var_ptr;

  GRBConstr *constr_ptr = d_model.getConstrs();
  for (size_t con = 0; con != d_m2; ++con)
    d_constrs.push_back(constr_ptr[con]);
  delete[] constr_ptr;
}