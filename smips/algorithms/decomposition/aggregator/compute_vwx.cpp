#include "aggregator.h"

vector<double> Aggregator::compute_vwx(double *x)
{
  vector<double> vx(d_probs.size());

  GRBConstr *cons = d_vw.getConstrs();  // move this in the loop for parallelization (use firstprivate(d_vw))
  GRBVar *vars = d_vw.getVars();

  for (size_t s = 0; s != d_probs.size(); ++s)
  {
    vector<double> rhs(d_problem.d_omega[s]);
    vector<vector<double>> &tech = d_problem.d_fix_tech ? d_problem.d_Tmat : d_problem.d_T_omega[s];
    for (size_t con = 0; con != rhs.size(); ++con)
      rhs[con] -= inner_product(tech[con].begin(), tech[con].end(), x, 0.0);
    d_vw.set(GRB_DoubleAttr_RHS, cons, rhs.data(), rhs.size());

    if (not d_fix_rec)    // update W and q
    {
      d_vw.set(GRB_DoubleAttr_Obj, vars, d_problem.d_q_omega[s].data(), d_problem.d_q_omega[s].size());
      vector<vector<double>> &Wmat = d_problem.d_W_omega[s];

      for (size_t row = 0; row != Wmat.size(); ++row)
      {
        vector<GRBConstr> constrs(Wmat[row].size(), cons[row]);
        d_vw.chgCoeffs(constrs.data(), vars, Wmat[row].data(), Wmat[row].size());
      }
    }
    d_vw.optimize();
    vx[s] = d_vw.get(GRB_DoubleAttr_ObjVal);
  }

  delete[] cons;
  delete[] vars;

  return vx;
}