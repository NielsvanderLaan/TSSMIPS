#include "sub.h"

void Sub::update(double *rhs, size_t s)
{
  d_model.set(GRB_DoubleAttr_RHS, d_constrs.data(), rhs, d_constrs.size());

  if (not d_problem.d_fix_rec)
  {
    d_model.set(GRB_DoubleAttr_Obj, d_vars.data(), d_problem.d_q_omega[s].data(), d_vars.size());



    for (size_t con = 0; con != d_constrs.size(); ++con)
    {
      vector<GRBConstr> con_vec(d_n2, d_constrs[con]);
      d_model.chgCoeffs(con_vec.data(), d_vars.data(), d_problem.d_W_omega[s][con].data(), d_vars.size());
    }


  }


}