#include "sub.h"


Sub::Sub(GRBEnv &env, Problem &problem, size_t s)
:
d_model(env),
d_problem(problem),
d_s(s)
{
  size_t n2 = problem.d_n2;
  size_t m2 = problem.d_m2;
  
  size_t ss_leq, ss_geq;
  ss_leq = problem.d_ss_leq; ss_geq = problem.d_ss_geq;
  
  d_n2 = n2;
  d_m2 = m2;

  vector<char> vTypes(n2, GRB_CONTINUOUS);
  vector<double> &q = problem.d_fix_rec ? problem.d_q : problem.d_q_omega[s];
  GRBVar* var_ptr = d_model.addVars(problem.d_l2.data(), problem.d_u2.data(), q.data(), vTypes.data(), NULL, n2);
  d_vars = vector<GRBVar> (var_ptr, var_ptr + n2);
  delete[] var_ptr;

      // constraint senses
  vector<char> senses(m2, GRB_EQUAL);
  fill(senses.begin(),          senses.begin() + ss_leq,          GRB_LESS_EQUAL);
  fill(senses.begin() + ss_leq, senses.begin() + ss_leq + ss_geq, GRB_GREATER_EQUAL);

  vector<double> rhs(m2, 0.0);

      // constraint lhs
  vector<vector<double>> &rm = problem.d_fix_rec ? problem.d_Wmat : problem.d_W_omega[s];
  GRBLinExpr Wy[m2];
  for (size_t conIdx = 0; conIdx != m2; ++conIdx)
  {
    double *row = rm[conIdx].data();
    Wy[conIdx].addTerms(row, d_vars.data(), n2);
  }
      // add constraints
  GRBConstr *constr_ptr = d_model.addConstrs(Wy, senses.data(), rhs.data(), NULL, m2);
  d_constrs = vector<GRBConstr> (constr_ptr, constr_ptr + m2);
  delete[] constr_ptr;

  d_model.update();
}