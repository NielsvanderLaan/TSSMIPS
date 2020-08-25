#include "lagrangian.h"

Lagrangian::Lagrangian(GRBEnv &env, Problem &problem)
:
d_model(env),
d_problem(problem)
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
  char zTypes[n1];
  fill_n(zTypes, p1, GRB_INTEGER);    
  fill_n(zTypes + p1, n1 - p1, GRB_CONTINUOUS);
  d_z_vars = d_model.addVars(problem.d_l1.data(), problem.d_u1.data(), NULL, zTypes, NULL, n1);   // cost coeffs set by update()

  vector<GRBLinExpr> lhs(problem.d_m1);
  for (size_t con = 0; con != problem.d_m1; ++con)
    lhs[con].addTerms(problem.d_Amat[con].data(), d_z_vars, n1);

  char senses1[lhs.size()];
  fill(senses1,                   senses1 + problem.d_fs_leq,          GRB_LESS_EQUAL);
  fill(senses1 + problem.d_fs_leq,          senses1 + problem.d_fs_leq + problem.d_fs_geq, GRB_GREATER_EQUAL);
  fill(senses1 + problem.d_fs_leq + problem.d_fs_geq, senses1 + lhs.size(),              GRB_EQUAL);
  delete[] d_model.addConstrs(lhs.data(), senses1, problem.d_b.data(), NULL, lhs.size());



  // adding second-stage variables (y)
       // variable types    
  char yTypes[n2];
  fill_n(yTypes, p2, GRB_INTEGER);    
  fill_n(yTypes + p2, n2 - p2, GRB_CONTINUOUS);
      // cost vector
  double *q = problem.d_q.data();        // transform cost vector and omega to c-style array 
      // add variables
  d_y_vars = d_model.addVars(problem.d_l2.data(), problem.d_u2.data(), q, yTypes, NULL, n2);

      // constraint senses
  char senses[m2];
  fill(senses,                   senses + ss_leq,          GRB_LESS_EQUAL);
  fill(senses + ss_leq,          senses + ss_leq + ss_geq, GRB_GREATER_EQUAL);
  fill(senses + ss_leq + ss_geq, senses + m2,              GRB_EQUAL);

      // constraint rhs
  double rhs[m2];
  fill(rhs, rhs + m2, 0.0);    

      // constraint lhs
  vector<vector<double>> &Wmat = problem.d_Wmat;
  vector<vector<double>> &Tmat = problem.d_Tmat;
  GRBLinExpr TxWy[m2];
  for (size_t conIdx = 0; conIdx != m2; ++conIdx)
  {
    TxWy[conIdx].addTerms(Tmat[conIdx].data(), d_z_vars, n1);
    TxWy[conIdx].addTerms(Wmat[conIdx].data(), d_y_vars, n2);
  }
      // add constraints
  GRBConstr *cons = d_model.addConstrs(TxWy, senses, rhs, NULL, m2);
  d_constrs = vector<GRBConstr> (cons, cons + m2);
  delete[] cons;
  d_model.update();
}