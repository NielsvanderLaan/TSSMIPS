#include "lagrangian.h"


Lagrangian::Lagrangian(GRBEnv &env, Problem &problem, size_t s)
:
d_model(env),
d_problem(problem),
d_rcut(false)
{
  size_t n1 = problem.d_n1;
  size_t p1 = problem.d_p1;
  size_t n2 = problem.d_n2;
  size_t p2 = problem.d_p2;
  size_t m2 = problem.d_m2;
  
  size_t ss_leq, ss_geq;
  ss_leq = problem.d_ss_leq; ss_geq = problem.d_ss_geq;
  
  d_n1 = n1;
  d_n2 = n2;
  d_m2 = m2;
  
  // adding first-stage variables (z)
  vector<char> zTypes(n1, GRB_CONTINUOUS);
  fill_n(zTypes.begin(), p1, GRB_INTEGER);
  GRBVar *zvars = d_model.addVars(problem.d_l1.data(), problem.d_u1.data(), NULL, zTypes.data(), NULL, zTypes.size());
  d_z_vars = vector<GRBVar>(zvars, zvars + n1);
  delete[] zvars;

  vector<GRBLinExpr> lhs(problem.d_m1);
  for (size_t con = 0; con != problem.d_m1; ++con)
    lhs[con].addTerms(problem.d_Amat[con].data(), d_z_vars.data(), n1);

  char senses1[lhs.size()];
  fill(senses1,                   senses1 + problem.d_fs_leq,          GRB_LESS_EQUAL);
  fill(senses1 + problem.d_fs_leq,          senses1 + problem.d_fs_leq + problem.d_fs_geq, GRB_GREATER_EQUAL);
  fill(senses1 + problem.d_fs_leq + problem.d_fs_geq, senses1 + lhs.size(),              GRB_EQUAL);
  delete[] d_model.addConstrs(lhs.data(), senses1, problem.d_b.data(), NULL, lhs.size());



  // adding second-stage variables (y)
       // variable types    
  vector<char> yTypes(n2, GRB_CONTINUOUS);
  fill_n(yTypes.begin(), p2, GRB_INTEGER);
  vector<double> &q = problem.d_fix_rec ? problem.d_q : problem.d_q_omega[s];
      // add variables
  GRBVar *yvars = d_model.addVars(problem.d_l2.data(), problem.d_u2.data(), q.data(), yTypes.data(), NULL, q.size());
  d_y_vars = vector<GRBVar>(yvars, yvars + n2);
  delete[] yvars;

      // constraint senses
  vector<char> senses(m2, GRB_EQUAL);
  fill_n(senses.begin(),          ss_leq, GRB_LESS_EQUAL);
  fill_n(senses.begin() + ss_leq, ss_geq, GRB_GREATER_EQUAL);

      // constraint rhs
  vector<double> &rhs = d_problem.d_omega[s];

      // constraint lhs
  vector<vector<double>> &Wmat = problem.d_fix_rec ? problem.d_Wmat : problem.d_W_omega[s];
  vector<vector<double>> &Tmat = problem.d_Tmat;
  GRBLinExpr TxWy[m2];
  for (size_t conIdx = 0; conIdx != m2; ++conIdx)
  {
    TxWy[conIdx].addTerms(Tmat[conIdx].data(), d_z_vars.data(), n1);
    TxWy[conIdx].addTerms(Wmat[conIdx].data(), d_y_vars.data(), n2);
  }
      // add constraints
  GRBConstr *cons = d_model.addConstrs(TxWy, senses.data(), rhs.data(), NULL, rhs.size());
  d_constrs = vector<GRBConstr> (cons, cons + m2);
  delete[] cons;

  d_theta = d_model.addVar(problem.d_L, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

  d_model.update();
}