#include "lagrangian.h"

void Lagrangian::update(size_t s, vector<double> &pi)
{
  d_model->set(GRB_DoubleAttr_RHS, d_constrs.data(), d_problem.d_omega[s].data(), d_constrs.size());

  if (not d_problem.d_fix_rec)
  {
    d_model->set(GRB_DoubleAttr_Obj, d_y_vars.data(), d_problem.d_q_omega[s].data(), d_n2);

    for (size_t con = 0; con != d_m2; ++con)
    {
      vector<GRBConstr> con_vec(d_n2, d_constrs[con]);
      d_model->chgCoeffs(con_vec.data(), d_y_vars.data(), d_problem.d_W_omega[s][con].data(), d_n2);
    }

  }

  update_pi(pi);
  d_model->update();
}